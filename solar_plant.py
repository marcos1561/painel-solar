import numpy as np
from numpy import cos, sin, arcsin, pi
import matplotlib.pyplot as plt
from matplotlib import cm
from typing import NamedTuple

from math_utils import rotation_matrix, numeric_integral

def solar_irradiance(cos_zenith_angle: np.ndarray):
    '''
    Irradiação solar em W/m^2 em função do ângulo zenital do Sol.

    Parâmetros:
    -----------
    cos_zenith_angle:
        cosseno do ângulo zenital do sol.
    '''
    irradiance = np.zeros(cos_zenith_angle.size)
    for i in range(cos_zenith_angle.size):
        if cos_zenith_angle[i] > 0.0001:
            am = 1/cos_zenith_angle[i]
            irradiance[i] = 1353 * 0.7**(am**0.678)
    
    return irradiance

class Panels(NamedTuple):
    efficiency: float
    num: int
    unit_area: float

class Plant():
    def __init__(self, a, b, w_r: float, w_t: float, e: float, 
        phi: float, panels_cfg: Panels) -> None:
        '''
        Parâmetros:
        -----------
        a e b: vetor 3-d (np.ndarray)
            Vértices de um quadrado de área 1 m^2, que possui a mesma orientação das
            placas solares. 
            
            Dado um sistema de coordenadas centrado no observador, considerando que
            um dos vértices desse quadrado esteja na origem desse sistema, sendo A
            esse vértice, a e b são vetores cuja cauda está em A e a ponta 
            está em um vértice adjacente à A.
            
            É calculado a energia incidente nesse
            quadrado, depois apenas é multiplicado o valor obtido pela área das placas.

            OBS: 
                -> O produto vetorial a X b deve dar o vetor que indica a orientação 
                da superfície da placa solar.

                -> a e b estão em metros.
        
        w_t:
            Velocidade angular de translação da terra em graus por dia.

        w_r:
            Velocidade angular de rotação da terra em graus por dia.

        e:
            Obliquidade da eclíptica em graus.

        phi:
            Latitude do observador em graus.
        '''
        self.a = a
        self.b = b
        self.w_t = w_t
        self.w_r = w_r
        self.e = e
        self.phi = phi
        self.w_s = solar_irradiance
        self.panels_cfg = panels_cfg

    def local_solar_pos(self, t: np.ndarray) -> np.ndarray:
        '''
        Dado os instantes de tempo em t (dia), calcula as respectivas posições do Sol
        no sistema de coordenadas do observador.

        Parâmetros:
            t:
                Array com os instates de tempo para calcular a posição do Sol no local em dias.

        Retorno:
            s_local: np.ndarray
                Array com shape=(3 , len(t)). A i-ésima coluna é a posição do Sol
                no local, no i-ésimo instante de tempo em "t".
        '''
        w_t = self.w_t
        w_r = self.w_r

        theta = t * w_t

        x_ecliptic = np.cos(theta)
        y_ecliptic = -np.sin(theta)
        z_ecliptic = np.zeros(t.size)

        s_ecliptic = np.array([x_ecliptic, y_ecliptic, z_ecliptic])

        s_equator = rotation_matrix([("y", self.e)]).dot(s_ecliptic)

        s_equator_rotated = []
        for idx, t_i in enumerate(t):
            alpha = w_r * t_i
            rot = rotation_matrix([("z", alpha)])
            s_i_local = rot.dot(s_equator[:, idx])
            s_equator_rotated.append(s_i_local)
        
        s_equator_rotated = np.array(s_equator_rotated).transpose()
        s_local = rotation_matrix([("x", np.pi/2 + self.phi)]).dot(s_equator_rotated)

        return s_local

    def proj_area(self, s_t: np.ndarray, verbose=False):
        s_t = s_t.transpose()
        
        s_dot_a = s_t.dot(self.a)
        proj_a = self.a - (s_dot_a * s_t.transpose()).transpose()
        
        s_dot_b = s_t.dot(self.b)
        proj_b = self.b - (s_dot_b * s_t.transpose()).transpose()

        cross_product = np.cross(proj_a, proj_b)
        proj_area_values = np.sqrt(cross_product[:, 0]**2 + cross_product[:, 1]**2 + cross_product[:,2]**2)

        proj_s_placa_scalar = s_t.dot(np.cross(self.a, self.b))

        below_placa = proj_s_placa_scalar < 0
        below_horizon = s_t[:,2] < 0

        proj_area_values[below_placa] = 0
        proj_area_values[below_horizon] = 0

        if verbose:
            '''
                Plota os vetores a e b (vértices da placa), suas projeções no espaço ortogonal a s
                (vetor que aponta para o Sol) e o vetor s.
            '''
            fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
            
            t_i = int(s_t.shape[0]/2)
            ax.quiver(0,0,0,*self.a, color="black", label="a")
            ax.quiver(0,0,0,*self.b, color="black", label="b")

            ax.quiver(0,0,0,*proj_a[t_i, :], color="green", label="proj_a")
            ax.quiver(0,0,0,*proj_b[t_i, :], color="red", label="proj_b")

            ax.quiver(0,0,0,*s_t[t_i, :], color="yellow", label="s")

            plt.show()

        return proj_area_values

    def proj_area_old(self, t, verbose=False):
        s_t = self.local_solar_pos(t)
        proj_a_to_s = np.ones(s_t.shape)
        proj_b_to_s = np.ones(s_t.shape)

        s_t_norm_sqr = s_t[0,:]**2 + s_t[1,:]**2 + s_t[2,:]**2
        proj_scalar_a = s_t.transpose().dot(self.a)/s_t_norm_sqr
        proj_scalar_b = s_t.transpose().dot(self.b)/s_t_norm_sqr

        proj_a_to_s[0,:] = proj_scalar_a * s_t[0,:]
        proj_a_to_s[1,:] = proj_scalar_a * s_t[1,:]
        proj_a_to_s[2,:] = proj_scalar_a * s_t[2,:]
        
        proj_b_to_s[0,:] = proj_scalar_b * s_t[0,:]
        proj_b_to_s[1,:] = proj_scalar_b * s_t[1,:]
        proj_b_to_s[2,:] = proj_scalar_b * s_t[2,:]

        proj_a = np.ones(s_t.shape)
        proj_b = np.ones(s_t.shape)

        proj_a[0,:] = self.a[0] - proj_a_to_s[0,:]
        proj_a[1,:] = self.a[1] - proj_a_to_s[1,:]
        proj_a[2,:] = self.a[2] - proj_a_to_s[2,:]
        
        proj_b[0,:] = self.b[0] - proj_b_to_s[0,:]
        proj_b[1,:] = self.b[1] - proj_b_to_s[1,:]
        proj_b[2,:] = self.b[2] - proj_b_to_s[2,:]

        cross_product = (np.cross(proj_a.transpose(), proj_b.transpose())).transpose()
        proj_area_values = np.sqrt(cross_product[0,:]**2 + cross_product[1,:]**2 + cross_product[2,:]**2)

        night = s_t[2,:] < 0
        
        placa_surface_orientation = np.cross(self.a, self.b)
        sun_down = s_t.transpose().dot(placa_surface_orientation) < 0

        without_sun = np.logical_or(night, sun_down)
        proj_area_values[without_sun] = 0

        if verbose:
            '''
                Plota os vetores a e b (vértices da placa), suas projeções no espaço ortogonal a s
                (vetor que aponta para o Sol) e o vetor s.
            '''
            fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
            
            t_i = int(t.size/2)
            ax.quiver(0,0,0,*self.a, color="black", label="a")
            ax.quiver(0,0,0,*self.b, color="black", label="b")

            ax.quiver(0,0,0,*proj_a[:,t_i], color="green", label="proj_a")
            ax.quiver(0,0,0,*proj_b[:,t_i], color="red", label="proj_b")

            ax.quiver(0,0,0,*s_t[:,t_i], color="yellow", label="s")

            plt.show()

        return proj_area_values

    def calculate_energy(self, t_i, t_f, delta_t, num_max_p=7000, progress=False):
        '''
        Calcula a energia gerada pela planta deste o instante `t_i` até `t_f` em KWh. Os instantes devem
        ser dados em dias, considerando que `t=0` é o Solstício de inverno no hemisfério Sul (máxima declinação) 
        com o ângulo horário do Sol igual a -90°, sendo que ângulo horário é 0° quando o sol está no meridiano local.

        Parâmetros:
        -----------
            t_i: float
                Tempo inicial da integração em dias.
            
            t_f: float
                Tempo final da integração em dias.
            
            delta_t: float
                Largura dos retângulos da integração em dias.

            num_max_p: int
                Número máximo de pontos utilizado em cada intervalo de tempo
                por passo da integração.
        
        Retorno:
            total_energy: float
                Energia gerada em KWh.
        '''
        t_range = np.arange(t_i, t_f+delta_t/2, delta_t)
        num_points = t_range.size

        count = 0
        num_steps = int(np.ceil((t_f-t_i)/delta_t / num_max_p))
        prog_freq = int(np.ceil(num_steps/15))
        
        total_integral = 0
        first_t_idx = 0
        while first_t_idx < num_points:
            last_t_idx = first_t_idx+num_max_p
            if last_t_idx > num_points:
                last_t_idx = num_points
            
            if progress:
                if count % prog_freq == 0:
                    print(f"Progresso: {count/num_steps*100:.2f} %\nComputando energia no intervalo ({t_range[first_t_idx]:.2f}, {t_range[last_t_idx-1]:.2f})\n")
            
            sub_t_range = t_range[first_t_idx: last_t_idx]

            s_t = self.local_solar_pos(sub_t_range)
            cos_zenith_angle = s_t[2]
            proj_area_values = self.proj_area(s_t)
            total_integral += numeric_integral(proj_area_values * self.w_s(cos_zenith_angle), delta_t)

            first_t_idx += num_max_p
            count += 1

        total_energy = total_integral * self.panels_cfg.num * self.panels_cfg.efficiency * self.panels_cfg.unit_area * 0.024 # Unidade está em Kwh

        return total_energy

def debug_visualization(planta: Plant, t_range):
    solar_pos = planta.local_solar_pos(t_range)
    area_values = planta.proj_area(solar_pos)
    # pot_values = area_values * planta.w_s * (planta.placas_cfg["total_area"]/planta.placa_model_area) * planta.placas_cfg["eff"] / 1000 # kW

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    ax.scatter(0, 0, 0, color="red")
    ax.plot([0, 1.2], [0, 0], [0, 0], color="red")
    ax.plot([0, 0], [0, 1.2], [0, 0], color="green")
    ax.plot([0, 0], [0, 0], [0, 1.2], color="black")

    ax.scatter(*planta.a, color="orange")
    ax.scatter(*planta.b, color="orange")
    ax.plot(*zip([0, 0, 0], list(planta.a)), linestyle='dashed', color="black", linewidth=2)
    ax.plot(*zip([0, 0, 0], list(planta.b)), linestyle='dashed', color="black", linewidth=2)

    theta_range = np.linspace(0, np.deg2rad(360), 100)
    ax.plot(np.cos(theta_range), np.sin(theta_range), color="black")

    ax.scatter(*solar_pos, c=np.linspace(0, 100, t_range.size), cmap="viridis")

    fig, ax = plt.subplots()
    ax.set_xlabel("Tempo (dias)")
    ax.set_ylabel("Potência (kW)")

    ax.plot(t_range, planta.panels_cfg.unit_area * area_values * planta.w_s(solar_pos[2]) * planta.panels_cfg.num * planta.panels_cfg.efficiency/1000)


    plt.show()
