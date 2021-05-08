from numpy import cos, sin, pi
import matplotlib.pyplot as plt
import functions as func

np = func.np

'''
Constantes
'''
phi = -28.93 / 180*pi # Latitude do local (negativo se for latitude Sul).
e = 23.4 / 180*pi # Obliquidade da eclíptica. 
w_s = 850 # Constante solar na superfície da Terra em W/s.
v_sa = 2*pi/365.25 # Velocidade angular do Sol em seu caminho na eclíptica em rad/dia.
v_sd = 2*pi # Velocidade angular do Sol na coordenada ângulo horário em rad/dia.


'''
Posição dos dois vértices, adjacentes ao vértice na origem, do quadrilátero.
'''
inclination = 27.5
p_1 = np.array([1,0,0])
p_2 = np.array([0, cos(inclination/180*pi), -sin(inclination/180*pi)])


'''
Limites de integração da integral da área projetada em função do tempo. A integração é feita numéricamente divindo a área em retângulos 
'''
delta_t = 0.1 # largura dos retângulos. 
t_x_max = 7000 # número máximo de retângulos que são calculados por vez.
t_i = 0 # tempo inicial (limite inferior)
t_f = 365.25  # tempo final (limite superior)


solar_energy = func.LigthArea(phi, e, w_s, v_sa, v_sd, p_1, p_2)

solar_energy.time_interval_steps(delta_t, t_x_max, t_i, t_f) # Definindo o intervalo de integração.
total_energy = solar_energy.calculate_energy() # Calculando a energia incidente
print(total_energy) 

# solar_energy.range_orientation([-4, 4], [27-1, 27+6], 1)
