from math import cos, sin, pi, radians
import datetime
import numpy as np

import solar_plant

'''
Constantes
'''
phi = radians(-38.93) # Latitude do local (negativo se for latitude Sul).
e = radians(23.4) # Obliquidade da eclíptica. 
# w_s = 850 # Constante solar na superfície da Terra em W.
w_t = radians(360)/365.25 # Velocidade angular do Sol em seu caminho na eclíptica em rad/dia.
w_r = radians(360) # Velocidade angular de rotação da Terra em rad/dia.
placas_cfg = solar_plant.Panels(
    efficiency=0.17,
    num=18,
    unit_area=2
)

'''
Posição dos dois vértices do quadrado, adjacentes ao vértice naa origem.
Esse quadrado deve ter área 1 e possuir a mesma orientação das placas solares.
'''
inclination = 20
# inclination = 37
p_1 = np.array([1,0,0])
p_2 = np.array([0, cos(inclination/180*pi), -sin(inclination/180*pi)])

'''
Intervalo de integração
'''
t_i, t_f, dt = 0, 360.25, 1/24/60*10
# t_i = (datetime.date(2022, 1, 20) - datetime.date(2022, 6, 21)).days 
# t_f = t_i + 1


planta = solar_plant.Plant(p_1, p_2, w_r, w_t, e, phi, placas_cfg)
# solar_plant.debug_visualization(planta, np.arange(t_i, t_f, dt))

energy = planta.calculate_energy(t_i, t_f, dt, 7000, progress=True)
print("Energia gerada:", round(energy,2), "KWh")
