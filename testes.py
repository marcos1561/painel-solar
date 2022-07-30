# import functions as utils
# import numpy as np
# import matplotlib.pyplot as plt


# phi = -30
# e = 23
# w_r = 360*0.9972
# w_t = 360/365.25


# def solar_path(t):
#     theta = np.deg2rad(w_t * t)
#     x = -np.sin(theta)
#     y = np.cos(theta)
#     z = np.zeros(t.size)

#     points = np.array([x, y, z]).transpose()

#     points_local = []
#     for idx, t_i in enumerate(t):
#         alpha = w_r * t_i
#         rot = utils.rotation_matrix([("x", np.deg2rad(e)), ("z", np.deg2rad(-alpha-90)), ("x", np.deg2rad(phi - 90))])

#         points_local.append(rot.dot(points[idx,:]))

#     points_local = np.array(points_local).transpose()
#     return points_local

# # x, y, z = line([0,0,1], [1,0,0], np.linspace(0, 10, 20))
# # length = np.sqrt(x**2 + y**2 + z**2)
# # x, y, z = x/length, y/length, z/length

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# ax.set_xlabel("x")
# ax.set_ylabel("y")
# ax.set_zlabel("z")

# ax.scatter(0, 0, 0, color="red")
# ax.plot([0, 1.2], [0, 0], [0, 0], color="red")
# ax.plot([0, 0], [0, 1.2], [0, 0], color="green")
# ax.plot([0, 0], [0, 0], [0, 1.2], color="black")

# theta_range = np.linspace(0, np.deg2rad(360), 100)
# ax.plot(np.cos(theta_range), np.sin(theta_range) , color="black")

# dia = 0
# t_range = np.linspace(dia, dia+1, 100)
# points = solar_path(t_range)

# ax.scatter(*points, c=np.linspace(0, 100, t_range.size), cmap="viridis")

# plt.show()