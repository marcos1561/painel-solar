from numpy import cos, sin, arcsin, pi
import numpy as np
import matplotlib.pyplot as plt
from numpy.core.arrayprint import printoptions
from numpy.core.fromnumeric import shape

phi = -28.93 / 180*pi
e = 23.4 / 180*pi
w_s = 850
v_sa = 2*pi/365.25
v_sd = 2*pi

p_1 = np.array([1,0,0])
p_2 = np.array([0, cos(30/180*pi), -0.5])

M_AB = np.array([[1,0,0], [0,sin(phi), cos(phi)], [0 , -cos(phi), sin(phi)]])


def d_t(t):
    return arcsin(-cos(v_sa*t)*sin(e))

def h_t(t):
    return -pi/2 + t*v_sd

def area_t(t):
    s_tb = np.array([-cos(d_t(t))*sin(h_t(t)), -cos(d_t(t))*cos(h_t(t)) , sin(d_t(t))])
    s_t = M_AB.dot(s_tb)

    n_t = t.size

    v_1 = np.array([np.ones(n_t), np.zeros(n_t), -s_t[0]/s_t[2]])
    
    v_norms = []
    for i in range(n_t):
        v_norms.append(1/np.linalg.norm(v_1[:,i]))

    v_norms = np.diag(v_norms)
    v_1 = v_1.dot(v_norms)

    v_2 = np.empty(shape=(3, n_t))
    for i in range(n_t):
        v_2[:,i] = np.cross(s_t[:,i], v_1[:,i])


    A = p_1.dot(v_1) * p_2.dot(v_2) - p_2.dot(v_1) * p_1.dot(v_2)
    return A, s_t[2]

delta_t = 0.01
t_x_max = 7000
t_i = 0
t_f = 365.25

size_step = (delta_t * t_x_max)
steps_t = np.floor( (t_f - t_i) / size_step )

ang_1arr = np.arange(-5 , 5, 0.5)
ang_n1 = len(ang_1arr)
ang_2arr = np.arange(27-5, 27+5, 0.5)
ang_n2 = len(ang_2arr)
print(ang_n1*ang_n2)

ang_1arr, ang_2arr = np.meshgrid(ang_1arr, ang_2arr)

energy = []
for ang_1, ang_2 in zip(ang_1arr.reshape(ang_n1*ang_n2), ang_2arr.reshape(ang_n1*ang_n2)):
    ang_1, ang_2 = ang_1/180*pi, ang_2/180*pi

    M_pr = np.array([[cos(ang_1), -sin(ang_1), 0], [sin(ang_1), cos(ang_1), 0], [0,0,1]])
    p_1 = M_pr.dot(np.array([1,0,0]))
    p_2 = M_pr.dot(np.array([0, cos(ang_2), -sin(ang_2)]))

    int_area = 0
    t_x = [t_i]
    for _ in range(0, int(steps_t)):
        int_area_step = 0
        t_x = np.arange(t_x[-1], t_x[-1] + size_step, delta_t)
        
        # print(f"{t_x[0]} | {t_x[-1]} | {len(t_x)}")

        area_s_z = area_t(t_x-delta_t/2)
        for area_i, s_z in zip(area_s_z[0], area_s_z[1]):
            if s_z > 0:
                if area_i > 0:
                    int_area_step += delta_t * area_i
        
        int_area += int_area_step

    if (size_step*(steps_t) < t_f):
        int_area_step = 0
        t_x = np.arange(t_x[-1], t_f, delta_t)
        
        # print(f"{t_x[0]} | {t_x[-1]} | {len(t_x)}")
        
        area_s_z = area_t(t_x-delta_t/2)
        for area_i, s_z in zip(area_s_z[0], area_s_z[1]):
            if s_z > 0:
                if area_i > 0:
                    int_area_step += delta_t * area_i
        
        int_area += int_area_step

    total_energy = int_area * w_s * 24 * 60 * 60 * 1 / 3600000
    print(f"{ang_1/pi*180};{ang_2/pi*180} | {total_energy}")
    energy.append(total_energy)


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

energy = np.array(energy).reshape(ang_n2, ang_n1)
np.savetxt("ang_1.txt", ang_1arr)
np.savetxt("ang_2.txt", ang_2arr)
np.savetxt("energy.txt", energy)

print(ang_1arr.reshape(ang_n2, ang_n1))
print(ang_2arr.reshape(ang_n2, ang_n1))
print(energy)

if len(energy) > 1:
    ax.plot_surface(ang_1arr, ang_2arr, energy)

    plt.show()

# t_x = np.linspace(t_i, t_f, 7000)
# area_tr = area_t(t_x)[0]

# plt.plot(t_x, area_tr)
# plt.xlabel("Tempo (dias)")
# plt.ylabel("Área capturada ($m^2$)")
# plt.grid()

# plt.show()
