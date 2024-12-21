"""
Defines symmetry groups that are useful for a variety of lattices.
"""

using Rotations
# Define the axes
x, y, z = [1, 0, 0], [0, 1, 0], [0, 0, 1]

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
    rotate_180yzm,
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
    rotate_2pi3mxmymz,
]

# C4 pi/2 rotation
rotate_90X = AngleAxis(π / 2, x...)
rotate_90Y = AngleAxis(π / 2, y...)
rotate_90Z = AngleAxis(π / 2, z...)
rotate_270X = AngleAxis(3 * π / 2, x...)
rotate_270Y = AngleAxis(3 * π / 2, y...)
rotate_270Z = AngleAxis(3 * π / 2, z...)
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

pyrochlore_symmetries = [ident, c2... , c2p... , c3... , c4... , inv, s4... , s6... , sih... , sid...]
