import numpy as np

def is_inside_baffled_wall(x, y, z, center, wall_radius, wall_thickness, wall_height, baffle_length, baffle_width):
    cx, cy, cz = center
    
    # Check Z bounds
    if z < cz - wall_height / 2.0 or z > cz + wall_height / 2.0:
        return False
        
    # Check Cylindrical Wall
    r = np.sqrt((x - cx)**2 + (y - cy)**2)
    r_inner = wall_radius
    r_outer = wall_radius + wall_thickness
    
    if r >= r_inner and r <= r_outer:
        return True
        
    # Check Baffles
    # Baffles are rectangles.
    # Baffle 1 (Right, 0 deg)
    # Extends from x = cx + r_inner - baffle_length to cx + r_inner
    # y range: cy - baffle_width/2 to cy + baffle_width/2
    if (x >= cx + r_inner - baffle_length and x <= cx + r_inner and
        y >= cy - baffle_width/2.0 and y <= cy + baffle_width/2.0):
        return True
        
    # Baffle 2 (Top, 90 deg)
    # Extends from y = cy + r_inner - baffle_length to cy + r_inner
    # x range: cx - baffle_width/2 to cx + baffle_width/2
    if (y >= cy + r_inner - baffle_length and y <= cy + r_inner and
        x >= cx - baffle_width/2.0 and x <= cx + baffle_width/2.0):
        return True
        
    # Baffle 3 (Left, 180 deg)
    # Extends from x = cx - r_inner to cx - r_inner + baffle_length
    # y range: cy - baffle_width/2 to cy + baffle_width/2
    if (x >= cx - r_inner and x <= cx - r_inner + baffle_length and
        y >= cy - baffle_width/2.0 and y <= cy + baffle_width/2.0):
        return True

    # Baffle 4 (Bottom, 270 deg)
    # Extends from y = cy - r_inner to cy - r_inner + baffle_length
    # x range: cx - baffle_width/2 to cx + baffle_width/2
    if (y >= cy - r_inner and y <= cy - r_inner + baffle_length and
        x >= cx - baffle_width/2.0 and x <= cx + baffle_width/2.0):
        return True
        
    return False

def generate_baffled_wall_csv(filename, nx, ny, nz, prob_lo, prob_hi, center, wall_radius, wall_thickness, wall_height, baffle_length, baffle_width):
    dx = (prob_hi[0] - prob_lo[0]) / nx
    dy = (prob_hi[1] - prob_lo[1]) / ny
    dz = (prob_hi[2] - prob_lo[2]) / nz
    
    print(f"Generating {filename}...")
    print(f"Grid: {nx}x{ny}x{nz}")
    print(f"Domain: {prob_lo} to {prob_hi}")
    
    with open(filename, 'w') as f:
        f.write("X,Y,Z,tag\n")
        for k in range(nz):
            z = prob_lo[2] + (k + 0.5) * dz
            for j in range(ny):
                y = prob_lo[1] + (j + 0.5) * dy
                for i in range(nx):
                    x = prob_lo[0] + (i + 0.5) * dx
                    
                    val = 0 # Fluid
                    if is_inside_baffled_wall(x, y, z, center, wall_radius, wall_thickness, wall_height, baffle_length, baffle_width):
                        val = 1 # Solid
                    
                    f.write(f"{i},{j},{k},{val}\n")
    print("Done.")

if __name__ == "__main__":
    nx, ny, nz = 660, 660, 4
    prob_lo = np.array([0.0, 0.0, -2.0])
    prob_hi = np.array([660.0, 660.0, 2.0])
    
    center = (330.0, 330.0, 0.0)
    wall_radius = 300.0
    wall_thickness = 400.0 # Updated thickness
    wall_height = 4.0
    
    baffle_length = 50.0
    baffle_width = 10.0
    
    generate_baffled_wall_csv("baffled_wall.csv", nx, ny, nz, prob_lo, prob_hi, center, wall_radius, wall_thickness, wall_height, baffle_length, baffle_width)
