import sys
import time
from pathlib import Path

def import_required_modules():
    """
    Import required dependencies and fail fast with actionable guidance if they
    are missing. The caller can install them explicitly (e.g. via `uv pip`
    or their preferred environment manager).
    """

    try:
        import numpy as np  # type: ignore
    except ImportError as exc:  # pragma: no cover - defensive
        raise RuntimeError(
            "Missing dependency 'numpy'. Install it via your environment manager,"
            " e.g. `uv pip install numpy`."
        ) from exc

    try:
        import tifffile as tf  # type: ignore
    except ImportError as exc:  # pragma: no cover - defensive
        raise RuntimeError(
            "Missing dependency 'tifffile'. Install it via your environment manager,"
            " e.g. `uv pip install tifffile`."
        ) from exc

    try:
        import imagecodecs  # type: ignore
        compression_available = True
    except ImportError:
        compression_available = False

    return np, tf, compression_available


def load_dimensions(dimensions_path: str):
    """Read the domain dimensions from the helper file."""

    with open(dimensions_path, "r", encoding="utf-8") as file_obj:
        dims = file_obj.read().strip().split()
    if len(dims) != 3:
        raise ValueError(
            f"Unexpected contents in {dimensions_path!r}: expected three integers"
        )
    return tuple(int(value) for value in dims)


def find_binary_file(n_x: int, n_y: int, n_z: int) -> Path:
    """Locate the binary file, allowing for optional suffixes like _seedX."""

    results_dir = Path("results")
    pattern = f"microstructure_nX{n_x}_nY{n_y}_nZ{n_z}*.bin"
    matches = sorted(results_dir.glob(pattern))

    if not matches:
        raise FileNotFoundError(
            f"No binary file matching pattern {pattern!r} in {results_dir.resolve()}"
        )

    if len(matches) > 1:
        chosen = matches[0]
        print(
            "Found multiple matching binary files, using the first one:"
            f" {chosen.name}\nCandidates: {[path.name for path in matches]}"
        )
        return chosen

    return matches[0]


def main() -> int:
    try:
        np, tf, compression_available = import_required_modules()
    except RuntimeError as exc:  # pragma: no cover - defensive
        print(exc, file=sys.stderr)
        return 1

    print("Starting binary to TIFF conversion...")
    start_time = time.time()

    n_x, n_y, n_z = load_dimensions("results/dimensions.txt")
    print(f"Dimensions: {n_x} x {n_y} x {n_z}")

    print("Loading binary data...")
    binary_path = find_binary_file(n_x, n_y, n_z)
    image_array = np.fromfile(str(binary_path), dtype=np.uint16)
    image_array = image_array.reshape((n_z, n_y, n_x))

    print("Saving TIFF...")
    tiff_path = binary_path.with_suffix(".tif")
    compression_args = {"compression": "lzw"} if compression_available else {}

    try:
        tf.imwrite(str(tiff_path), image_array, **compression_args)
        if compression_available:
            print("Saved with LZW compression")
        else:
            print("Saved without compression")
    except Exception as exc:  # pragma: no cover - defensive
        print(f"Compression failed: {exc}")
        print("Retrying without compression...")
        tf.imwrite(str(tiff_path), image_array)
        print("Saved without compression")

    end_time = time.time()
    print(f"TIFF conversion completed in {end_time - start_time:.2f} seconds")
    print(f"TIFF file saved: {tiff_path}")
    return 0


if __name__ == "__main__":  # pragma: no cover - script entry point
    sys.exit(main())