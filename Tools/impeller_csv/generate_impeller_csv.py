import numpy as np

def point_in_polygon(x, y, poly):
    """
    Ray casting algorithm to check if point (x,y) is inside polygon.
    poly: list of (x,y) tuples
    """
    n = len(poly)
    inside = False
    p1x, p1y = poly[0]
    for i in range(n + 1):
        p2x, p2y = poly[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y
    return inside

def generate_impeller_csv(filename, nx, ny, nz, prob_lo, prob_hi, center, hub_radius, blade_outer_radius, blade_thickness, height, num_blades=8):
    dx = (prob_hi[0] - prob_lo[0]) / nx
    dy = (prob_hi[1] - prob_lo[1]) / ny
    dz = (prob_hi[2] - prob_lo[2]) / nz
    
    cx, cy, cz = center
    
    # Generate Polygon Profile (Same logic as STL generator)
    profile_points = [] 
    alpha = np.arcsin((blade_thickness / 2.0) / hub_radius)
    
    for k in range(num_blades):
        phi = 2 * np.pi * k / num_blades
        next_phi = 2 * np.pi * (k + 1) / num_blades
        
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        
        x_root = np.sqrt(hub_radius**2 - (blade_thickness/2.0)**2)
        
        p_root_right_local = np.array([x_root, -blade_thickness/2.0])
        p_tip_right_local  = np.array([blade_outer_radius, -blade_thickness/2.0])
        p_tip_left_local   = np.array([blade_outer_radius, blade_thickness/2.0])
        p_root_left_local  = np.array([x_root, blade_thickness/2.0])
        
        rot_mat = np.array([[cos_phi, -sin_phi], [sin_phi, cos_phi]])
        
        profile_points.append(rot_mat @ p_root_right_local)
        profile_points.append(rot_mat @ p_tip_right_local)
        profile_points.append(rot_mat @ p_tip_left_local)
        profile_points.append(rot_mat @ p_root_left_local)
        
        start_angle = phi + alpha
        end_angle = next_phi - alpha
        
        num_arc_points = 5
        for j in range(num_arc_points):
            theta = start_angle + (end_angle - start_angle) * (j + 1) / (num_arc_points + 1)
            px = hub_radius * np.cos(theta)
            py = hub_radius * np.sin(theta)
            profile_points.append(np.array([px, py]))

    # Shift profile to center
    poly = [(cx + p[0], cy + p[1]) for p in profile_points]
    
    z_min = cz - height / 2.0
    z_max = cz + height / 2.0

    print(f"Generating {filename}...")
    print(f"Grid: {nx}x{ny}x{nz}")
    print(f"Domain: {prob_lo} to {prob_hi}")
    print(f"Impeller: Center={center}")

    with open(filename, 'w') as f:
        f.write("X,Y,Z,tag\n")
        for k in range(nz):
            z = prob_lo[2] + (k + 0.5) * dz
            in_z = (z >= z_min) and (z <= z_max)
            
            for j in range(ny):
                y = prob_lo[1] + (j + 0.5) * dy
                for i in range(nx):
                    x = prob_lo[0] + (i + 0.5) * dx
                    
                    val = 0 # Fluid
                    if in_z:
                        if point_in_polygon(x, y, poly):
                            val = 1 # Solid
                    
                    f.write(f"{i},{j},{k},{val}\n")
    print("Done.")

if __name__ == "__main__":
    nx, ny, nz = 660, 660, 4
    prob_lo = np.array([0.0, 0.0, -2.0])
    prob_hi = np.array([660.0, 660.0, 2.0])
    
    center = (330.0, 330.0, 0.0)
    hub_radius = 10.0
    blade_outer_radius = 15.0
    blade_thickness = 2.0
    height = 10.0
    
    generate_impeller_csv("impeller.csv", nx, ny, nz, prob_lo, prob_hi, center, hub_radius, blade_outer_radius, blade_thickness, height)
