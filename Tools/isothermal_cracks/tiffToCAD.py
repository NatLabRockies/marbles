import sys
import time
import os
import traceback


def import_required_modules():
    """
    Import the Python dependencies required for TIFF to CAD conversion.

    The function raises a RuntimeError with installation guidance if a
    dependency is missing so the caller can install it explicitly (for
    example via `uv pip install <package>`).
    """

    try:
        import numpy as np  # type: ignore
    except ImportError as exc:  # pragma: no cover - defensive
        raise RuntimeError(
            "Missing dependency 'numpy'. Install it, e.g. `uv pip install numpy`."
        ) from exc

    try:
        import tifffile as tf  # type: ignore
    except ImportError as exc:  # pragma: no cover - defensive
        raise RuntimeError(
            "Missing dependency 'tifffile'. Install it, e.g. `uv pip install tifffile`."
        ) from exc

    try:
        from skimage import measure  # type: ignore
    except ImportError as exc:  # pragma: no cover - defensive
        raise RuntimeError(
            "Missing dependency 'scikit-image' (provides skimage.measure)."
            " Install it, e.g. `uv pip install scikit-image`."
        ) from exc

    try:
        import trimesh  # type: ignore
    except ImportError as exc:  # pragma: no cover - defensive
        raise RuntimeError(
            "Missing dependency 'trimesh'. Install it, e.g. `uv pip install trimesh`."
        ) from exc

    try:
        import FreeCAD  # type: ignore
        import Part  # type: ignore
        import Mesh  # type: ignore
        freecad_available = True
    except ImportError:
        FreeCAD = None
        Part = None
        Mesh = None
        freecad_available = False

    return {
        "np": np,
        "tf": tf,
        "measure": measure,
        "trimesh": trimesh,
        "freecad_available": freecad_available,
        "FreeCAD": FreeCAD,
        "Part": Part,
        "Mesh": Mesh,
    }


def find_tiff_file():
    """Find the TIFF file in results directory or current directory"""
    # Check current directory first
    current_dir = "."
    tiff_files = [f for f in os.listdir(current_dir) if f.endswith('.tif')]
    
    if tiff_files:
        tiff_file = tiff_files[0]
        if len(tiff_files) > 1:
            print(f"Multiple TIFF files found in current directory: {tiff_files}")
            print(f"Using: {tiff_file}")
        return tiff_file
    
    # Check results directory
    results_dir = "results"
    if os.path.exists(results_dir):
        tiff_files = [f for f in os.listdir(results_dir) if f.endswith('.tif')]
        if tiff_files:
            tiff_file = os.path.join(results_dir, tiff_files[0])
            if len(tiff_files) > 1:
                print(f"Multiple TIFF files found in results directory: {tiff_files}")
                print(f"Using: {tiff_files[0]}")
            return tiff_file
    
    raise FileNotFoundError("No TIFF files found in current directory or results/")

def generate_mesh_for_tag(
    volume_data,
    tag_value,
    np_module,
    measure_module,
    trimesh_module,
    smoothing=True,
):
    """Generate mesh for a specific tag value"""
    if tag_value == -1:
        # Special case for combined mesh
        binary_mask = (volume_data > 0).astype(np_module.uint8)
        print(f"Generating combined mesh for all solid materials...")
    else:
        print(f"Generating mesh for tag {tag_value}...")
        binary_mask = (volume_data == tag_value).astype(np_module.uint8)
    
    # Check if tag exists
    if np_module.sum(binary_mask) == 0:
        print(f"Warning: No voxels found for tag {tag_value if tag_value != -1 else 'combined'}")
        return None
    
    print(
        f"Found {np_module.sum(binary_mask)} voxels for tag "
        f"{tag_value if tag_value != -1 else 'combined'}"
    )
    
    # Generate mesh using marching cubes
    try:
        verts, faces, normals, values = measure_module.marching_cubes(
            binary_mask, 
            level=0.5,
            spacing=(1.0, 1.0, 1.0)
        )
        
        print(f"Generated mesh: {len(verts)} vertices, {len(faces)} faces")
        
        # Create trimesh object
        mesh = trimesh_module.Trimesh(vertices=verts, faces=faces, vertex_normals=normals)
        
        # Apply smoothing using correct trimesh API
        if smoothing and len(verts) > 0:
            try:
                # Use smoothing filter if available
                mesh = mesh.smoothed()
            except AttributeError:
                try:
                    # Alternative smoothing method
                    mesh = trimesh_module.smoothing.filter_laplacian(mesh, iterations=1)
                except (AttributeError, ImportError):
                    # Skip smoothing if not available
                    print("Smoothing not available - using original mesh")
        
        # Clean up mesh
        mesh.remove_degenerate_faces()
        mesh.remove_duplicate_faces()
        mesh.remove_unreferenced_vertices()
        
        if not mesh.is_empty:
            print(f"Final mesh: {len(mesh.vertices)} vertices, {len(mesh.faces)} faces")
            return mesh
        else:
            print(f"Warning: Empty mesh generated for tag {tag_value}")
            return None
            
    except Exception as e:
        print(f"Error generating mesh for tag {tag_value}: {e}")
        return None

def export_stl(mesh, filename, tag_value):
    """Export mesh to STL format"""
    try:
        if tag_value == -1:
            stl_filename = f"{filename}_combined.stl"
        else:
            stl_filename = f"{filename}_tag{tag_value}.stl"
        mesh.export(stl_filename)
        print(f"Exported STL: {stl_filename}")
        return stl_filename
    except Exception as e:
        print(f"Error exporting STL for tag {tag_value}: {e}")
        return None

def export_step_iges_freecad(
    mesh,
    filename,
    tag_value,
    freecad_modules,
    format='step',
):
    """Export mesh to STEP/IGES using FreeCAD"""
    FreeCAD = freecad_modules[0]
    Part = freecad_modules[1]
    Mesh_module = freecad_modules[2]

    if FreeCAD is None or Part is None or Mesh_module is None:
        print(f"FreeCAD not available - cannot export {format.upper()}")
        return None
    
    try:
        # Create FreeCAD document
        doc = FreeCAD.newDocument()
        
        # Convert trimesh to FreeCAD mesh
        freecad_mesh = Mesh_module.Mesh()
        
        # Add triangles to FreeCAD mesh
        for face in mesh.faces:
            triangle = [mesh.vertices[face[0]], mesh.vertices[face[1]], mesh.vertices[face[2]]]
            freecad_mesh.addFacet(triangle[0], triangle[1], triangle[2])
        
        # Create mesh object in document
        mesh_obj = doc.addObject("Mesh::Feature", f"Tag{tag_value}")
        mesh_obj.Mesh = freecad_mesh
        
        # Export based on format
        if format.lower() == 'step':
            if tag_value == -1:
                output_filename = f"{filename}_combined.step"
            else:
                output_filename = f"{filename}_tag{tag_value}.step"
            Part.export([mesh_obj], output_filename)
                
        elif format.lower() == 'iges':
            if tag_value == -1:
                output_filename = f"{filename}_combined.iges"
            else:
                output_filename = f"{filename}_tag{tag_value}.iges"
            Part.export([mesh_obj], output_filename)
        else:
            raise ValueError(f"Unsupported export format: {format}")
        
        # Close document
        FreeCAD.closeDocument(doc.Name)
        
        print(f"Exported {format.upper()}: {output_filename}")
        return output_filename
        
    except Exception as e:
        print(f"Error exporting {format.upper()} for tag {tag_value}: {e}")
        return None

def main() -> int:
    print("Starting TIFF to CAD conversion...")
    start_time = time.time()

    try:
        modules = import_required_modules()
    except RuntimeError as exc:  # pragma: no cover - defensive
        print(exc, file=sys.stderr)
        return 1

    np_module = modules["np"]
    tf_module = modules["tf"]
    measure_module = modules["measure"]
    trimesh_module = modules["trimesh"]
    freecad_available = modules["freecad_available"]
    freecad_modules = (
        modules["FreeCAD"],
        modules["Part"],
        modules["Mesh"],
    )

    if freecad_available:
        print("FreeCAD available for STEP/IGES export")
    else:
        print("FreeCAD not available - only STL export supported")

    try:
        tiff_file = find_tiff_file()
        print(f"Loading TIFF file: {tiff_file}")

        volume_data = tf_module.imread(tiff_file)
        print(f"TIFF dimensions: {volume_data.shape}")
        print(f"Data type: {volume_data.dtype}")

        unique_tags = np_module.unique(volume_data)
        unique_tags = unique_tags[unique_tags > 0]
        print(f"Found tags: {unique_tags}")

        if len(unique_tags) == 0:
            print("No non-zero tags found in TIFF file")
            return 0

        output_dir = "cad_exports"
        os.makedirs(output_dir, exist_ok=True)

        base_filename = os.path.splitext(os.path.basename(tiff_file))[0]

        exported_files = []

        max_tags_to_process = min(100, len(unique_tags))
        print(
            f"Processing first {max_tags_to_process} tags (out of {len(unique_tags)})"
            " for faster testing..."
        )

        for index, tag in enumerate(unique_tags[:max_tags_to_process]):
            print(f"\nProcessing tag {tag} ({index + 1}/{max_tags_to_process})...")

            mesh = generate_mesh_for_tag(
                volume_data,
                tag,
                np_module,
                measure_module,
                trimesh_module,
                smoothing=False,
            )

            if mesh is not None:
                output_base = os.path.join(output_dir, base_filename)

                stl_file = export_stl(mesh, output_base, tag)
                if stl_file:
                    exported_files.append(stl_file)

                if freecad_available:
                    step_file = export_step_iges_freecad(
                        mesh,
                        output_base,
                        tag,
                        freecad_modules,
                        "step",
                    )
                    if step_file:
                        exported_files.append(step_file)

                    iges_file = export_step_iges_freecad(
                        mesh,
                        output_base,
                        tag,
                        freecad_modules,
                        "iges",
                    )
                    if iges_file:
                        exported_files.append(iges_file)
                else:
                    print("For STEP/IGES export, install FreeCAD (e.g. via conda).")

        print(f"\nGenerating combined mesh...")
        combined_mesh = generate_mesh_for_tag(
            volume_data,
            -1,
            np_module,
            measure_module,
            trimesh_module,
            smoothing=False,
        )

        if combined_mesh is not None:
            output_base = os.path.join(output_dir, base_filename)

            stl_file = export_stl(combined_mesh, output_base, -1)
            if stl_file:
                exported_files.append(stl_file)

            if freecad_available:
                step_file = export_step_iges_freecad(
                    combined_mesh,
                    output_base,
                    -1,
                    freecad_modules,
                    "step",
                )
                if step_file:
                    exported_files.append(step_file)

                iges_file = export_step_iges_freecad(
                    combined_mesh,
                    output_base,
                    -1,
                    freecad_modules,
                    "iges",
                )
                if iges_file:
                    exported_files.append(iges_file)

        end_time = time.time()
        print(f"\nCAD conversion completed in {end_time - start_time:.2f} seconds")
        print(f"Exported {len(exported_files)} files:")
        for file_path in exported_files:
            if os.path.exists(file_path):
                file_size = os.path.getsize(file_path) / (1024 * 1024)
                print(f"  - {file_path} ({file_size:.1f} MB)")
            else:
                print(f"  - {file_path} (file not found)")

        print(f"\nFiles saved in: {output_dir}/")
        print("\nCAD files can be opened in:")
        print("- STL: Blender, MeshLab, FreeCAD, SolidWorks, Fusion 360")
        print("- STEP: FreeCAD, SolidWorks, Fusion 360, Onshape")
        print("- IGES: FreeCAD, SolidWorks, Fusion 360")

        if len(unique_tags) > max_tags_to_process:
            print(
                f"\nNote: Only processed {max_tags_to_process} out of {len(unique_tags)} tags."
            )
            print("To process all tags, change max_tags_to_process in the script.")

        return 0

    except Exception as exc:  # pragma: no cover - defensive
        print(f"Error during conversion: {exc}")
        traceback.print_exc()
        return 1

if __name__ == "__main__":  # pragma: no cover - script entry point
    sys.exit(main())