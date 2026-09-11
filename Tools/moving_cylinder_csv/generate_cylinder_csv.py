import numpy as np

def generate_cylinder_csv(filename, nx, ny, nz, prob_lo, prob_hi, cyl_center, cyl_radius):
    dx = (prob_hi[0] - prob_lo[0]) / nx
    dy = (prob_hi[1] - prob_lo[1]) / ny
    dz = (prob_hi[2] - prob_lo[2]) / nz
    
    print(f"Generating {filename}...")
    print(f"Grid: {nx}x{ny}x{nz}")
    print(f"Domain: {prob_lo} to {prob_hi}")
    print(f"Cylinder: Center={cyl_center}, Radius={cyl_radius}")

    # Open file
    with open(filename, 'w') as f:
        f.write("X,Y,Z,tag\n")
        # Loop order must match AMReX: k, j, i (Z, Y, X)
        for k in range(nz):
            z = prob_lo[2] + (k + 0.5) * dz
            for j in range(ny):
                y = prob_lo[1] + (j + 0.5) * dy
                for i in range(nx):
                    x = prob_lo[0] + (i + 0.5) * dx
                    
                    # Check if inside cylinder (Z-aligned)
                    dist_sq = (x - cyl_center[0])**2 + (y - cyl_center[1])**2
                    
                    if dist_sq <= cyl_radius**2:
                        val = 1 # Solid
                    else:
                        val = 0 # Fluid
                    
                    f.write(f"{i},{j},{k},{val}\n")
    print("Done.")

if __name__ == "__main__":
    # Cylinder Case
    nx, ny, nz = 660, 123, 4
    prob_lo = np.array([0.0, 0.0, -2.0])
    prob_hi = np.array([660.0, 123.0, 2.0])
    
    cyl_center = (600.0, 61.5, 0.0)
    cyl_radius = 15.0
    
    generate_cylinder_csv("cylinder.csv", nx, ny, nz, prob_lo, prob_hi, cyl_center, cyl_radius)
