module ex4

using NLCE
using Rotations

basis = [[0, 0, 0], [1 / 4, 1 / 4, 1 / 4]]

primitive_vec = [[0, 1 / 2, 1 / 2], [1 / 2, 0, 1 / 2], [1 / 2, 1 / 2, 0]]

neighborhood = [sqrt(3) / 4]

# Try for pyrochore instead
basis = [[0, 0, 0], [0, 1 / 4, 1 / 4], [1 / 4, 0, 1 / 4], [1 / 4, 1 / 4, 0]]

primitive_vec = [[0, 1 / 2, 1 / 2], [1 / 2, 0, 1 / 2], [1 / 2, 1 / 2, 0]]

neighborhood = [sqrt(2) / 4]
# Setting the maximum order
max_order = 4

x, y, z = [1, 0, 0], [0, 1, 0], [0, 0, 1]

# Seems to be okay with Identity, pi rotations about x, y, z, 
#[2916, 1957, 1957, 1957, 1, 1, 1, 1, 1, 1, 2916, 1957, 1957, 1957, 1957, 1957, 1957, 2916, 1, 1, 1, 1, 1, 1, 1, 1957, 1957, 1957, 1957, 1957, 1957, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1957, 1957, 1957, 2916, 2916, 2916]
# Identity
ident = [1 0 0; 0 1 0; 0 0 1]

# C2 - pi rotations about x, y and z axes
rotate_180X = AngleAxis(π, x...)
rotate_180Y = AngleAxis(π, y...)
rotate_180Z = AngleAxis(π, z...)
c2 = [rotate_180X, rotate_180Y, rotate_180Z]

# C2' - pi rotation about x + y, x + z, y + z, cubic face diagonals
rotate_180xyp = AngleAxis(π, (x + y)...)
rotate_180xzp = AngleAxis(π, (x + z)...)
rotate_180yzp = AngleAxis(π, (z + y)...)
rotate_180xym = AngleAxis(π, (x - y)...)
rotate_180xzm = AngleAxis(π, (x - z)...)
rotate_180yzm = AngleAxis(π, (y - z)...)
c2p = [
        rotate_180xyp,
        rotate_180xzp,
        rotate_180yzp,
        rotate_180xym,
        rotate_180xzm,
        rotate_180yzm
        ]

# C3 - 2pi/3 rotations about 1, 1, 1 cubic body diagonals
rotate_2pi3pxpypz = AngleAxis(2 * π / 3, (x + y + z)...)
rotate_2pi3pxpymz = AngleAxis(2 * π / 3, (x + y + -z)...)
rotate_2pi3pxmypz = AngleAxis(2 * π / 3, (x + -y + z)...)
rotate_2pi3mxpypz = AngleAxis(2 * π / 3, (-x + y + z)...)
rotate_2pi3mxmypz = AngleAxis(2 * π / 3, (-x + -y + z)...)
rotate_2pi3mxpymz = AngleAxis(2 * π / 3, (-x + y + -z)...)
rotate_2pi3pxmymz = AngleAxis(2 * π / 3, (x + -y + -z)...)
rotate_2pi3mxmymz = AngleAxis(2 * π / 3, (-x + -y + -z)...)
c3 = [
    rotate_2pi3pxpypz,
    rotate_2pi3pxpymz,
    rotate_2pi3pxmypz,
    rotate_2pi3mxpypz,
    rotate_2pi3mxmypz,
    rotate_2pi3mxpymz,
    rotate_2pi3pxmymz,
    rotate_2pi3mxmymz
    ]

# C4 pi/2 rotation
rotate_90X = AngleAxis(π/2, x...)
rotate_90Y = AngleAxis(π/2, y...)
rotate_90Z = AngleAxis(π/2, z...)
rotate_270X = AngleAxis(3 * π/2, x...)
rotate_270Y = AngleAxis(3 * π/2, y...)
rotate_270Z = AngleAxis(3 * π/2, z...)
c4 = [rotate_90X, rotate_90Y, rotate_90Z, rotate_270X, rotate_270Y, rotate_270Z]

# Inversion operator
inv = -[1 0 0; 0 1 0; 0 0 1]

# Inversion times group operators
# S4 pi/2 rotation about x, y, z, then inversion
s4 = [inv * x for x in c4]

# S6 pi/3 rotations about 1 1 1 cubic body diagonal, then reflection
s6 = [inv * x for x in c3]

# sigma h reflection through the planes normal to the C4 rotations (rotate by c2, then flip)
sih = [inv * x for x in c2]

# sigma d reflections through c2' rotations
sid = [inv * x for x in c2p]

group = [[ident], c2, c2p, c3, c4, [inv], s4, s6, sih, sid]
flat_group = []

for elem_arr in group
    for elem in elem_arr
        push!(flat_group, elem)
    end
end

#coordinates, colors, centers = NLCE.generate_coordinates(basis, primitive_vec, 2, repeat([1], length(basis)))
##println(coordinates)
#println(NLCE.find_permutations(coordinates, basis, flat_group[1:3]))

# Generating all the clusters using this information
nlce_clusters = coord_NLCE(flat_group, basis, primitive_vec, neighborhood, max_order)

# Writing all the files to the corresponding folder, creating the folder
# if it does not exist
filepath = "examples/outputs/ex-5/pyrochlore_nn_sym"
mkpath(filepath)
filename = filepath * "/pyrochlore_nn_sym"

# Write all the files in the default format
NLCE.write_to_file_coordinates(nlce_clusters, filename)

end
