using GTO, LinearAlgebra
c_1 = CartesianGaussian(0, 0, 0, 1, [0, 0, 0])
c_2 = CartesianGaussian(0, 0, 0, 1, [1, 0, 0])

pos_integral(c_1, c_1)
overlap_integral(c_2 * c_2)
pos_integral(c_2, c_2)
c_2 = CartesianGaussian(1, 2, 4, 1, [2, 0, 0]);
r2 = r2_integral(c_2, c_2) / overlap_integral(c_2 * c_2);
r = pos_integral(c_2, c_2) / overlap_integral(c_2 * c_2);
r2 - norm(r)^2

