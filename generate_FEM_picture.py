import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Create a triangular mesh
def create_triangular_mesh(nx=5, ny=5):
    """Create a simple triangular mesh on a rectangular domain"""
    x = np.linspace(0, 4, nx)
    y = np.linspace(0, 3, ny)
    X, Y = np.meshgrid(x, y)
    
    # Flatten for vertices
    vertices = np.column_stack([X.ravel(), Y.ravel()])
    
    # Create triangles
    triangles = []
    for i in range(ny - 1):
        for j in range(nx - 1):
            # Lower triangle
            v1 = i * nx + j
            v2 = i * nx + j + 1
            v3 = (i + 1) * nx + j
            triangles.append([v1, v2, v3])
            
            # Upper triangle
            v1 = i * nx + j + 1
            v2 = (i + 1) * nx + j + 1
            v3 = (i + 1) * nx + j
            triangles.append([v1, v2, v3])
    
    return vertices, np.array(triangles)

# Linear basis function
def basis_function(vertices, triangles, node_idx):
    """Compute the linear basis function value at each vertex"""
    z = np.zeros(len(vertices))
    z[node_idx] = 1.0
    return z

# Visualization
fig = plt.figure(figsize=(14, 10))
ax = fig.add_subplot(111, projection='3d')

# Create mesh
vertices, triangles = create_triangular_mesh(6, 5)

# Choose two nodes to highlight (e.g., nodes i and j)
node_i = 14  # Adjust based on your mesh
node_j = 15  # Adjust based on your mesh

# Compute basis functions
z_i = basis_function(vertices, triangles, node_i)
z_j = basis_function(vertices, triangles, node_j)

# Plot the mesh edges
for tri in triangles:
    pts = vertices[tri]
    # Create edges
    edges = [[pts[0], pts[1]], [pts[1], pts[2]], [pts[2], pts[0]]]
    for edge in edges:
        xs = [edge[0][0], edge[1][0]]
        ys = [edge[0][1], edge[1][1]]
        zs = [0, 0]
        ax.plot(xs, ys, zs, 'gray', alpha=0.3, linewidth=0.5)

# Plot triangular elements with basis function i
for tri in triangles:
    pts = vertices[tri]
    z_vals_i = z_i[tri]
    
    if np.any(z_vals_i > 0):  # Only plot if this triangle is in the support
        verts = [list(zip(pts[:, 0], pts[:, 1], z_vals_i))]
        poly = Poly3DCollection(verts, alpha=0.6, facecolor='blue', edgecolor='darkblue', linewidth=0.5)
        ax.add_collection3d(poly)

# Plot triangular elements with basis function j
for tri in triangles:
    pts = vertices[tri]
    z_vals_j = z_j[tri]
    
    if np.any(z_vals_j > 0):  # Only plot if this triangle is in the support
        verts = [list(zip(pts[:, 0], pts[:, 1], z_vals_j))]
        poly = Poly3DCollection(verts, alpha=0.6, facecolor='mediumpurple', edgecolor='purple', linewidth=0.5)
        ax.add_collection3d(poly)

# Mark the special nodes
ax.scatter([vertices[node_i, 0]], [vertices[node_i, 1]], [0], 
           color='black', s=100, zorder=5)
ax.scatter([vertices[node_j, 0]], [vertices[node_j, 1]], [0], 
           color='black', s=100, zorder=5)

# Draw vertical lines to the peaks
ax.plot([vertices[node_i, 0], vertices[node_i, 0]], 
        [vertices[node_i, 1], vertices[node_i, 1]], 
        [0, 1], 'k--', linewidth=1.5)
ax.plot([vertices[node_j, 0], vertices[node_j, 0]], 
        [vertices[node_j, 1], vertices[node_j, 1]], 
        [0, 1], 'k--', linewidth=1.5)

# Mark the peaks
ax.scatter([vertices[node_i, 0]], [vertices[node_i, 1]], [1], 
           color='darkblue', s=100, zorder=5)
ax.scatter([vertices[node_j, 0]], [vertices[node_j, 1]], [1], 
           color='purple', s=100, zorder=5)

# Add labels
ax.text(vertices[node_i, 0], vertices[node_i, 1], -0.15, 'Node i', 
        fontsize=12, ha='center')
ax.text(vertices[node_j, 0], vertices[node_j, 1], -0.15, 'Node j', 
        fontsize=12, ha='center')
ax.text(vertices[node_i, 0], vertices[node_i, 1], 1.15, r'$\psi_i$', 
        fontsize=14, ha='center')
ax.text(vertices[node_j, 0], vertices[node_j, 1], 1.15, r'$\psi_j$', 
        fontsize=14, ha='center')

"""
# Add title
ax.text2D(0.5, 0.95, 'Basis Function exampla', 
          transform=ax.transAxes, fontsize=14, ha='center')
ax.text2D(0.85, 0.05, 'Triangular elements', 
          transform=ax.transAxes, fontsize=10, ha='center')
"""
# Set view angle
ax.view_init(elev=22.5, azim=62.5)

# Set labels and limits
ax.set_xlabel('x', fontsize=10)
ax.set_ylabel('y', fontsize=10)
ax.set_zlabel('z', fontsize=10)
ax.set_zlim(-0.2, 1.3)

# Remove axis ticks for cleaner look
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

# Set background color
ax.xaxis.pane.fill = True
ax.yaxis.pane.fill = True
ax.zaxis.pane.fill = True
ax.xaxis.pane.set_facecolor('#E8F4F8')
ax.yaxis.pane.set_facecolor('#E8F4F8')
ax.zaxis.pane.set_facecolor('#E8F4F8')

plt.tight_layout()
plt.show()