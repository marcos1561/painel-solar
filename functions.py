from numpy import cos, sin, arcsin, pi
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

class LigthArea(object):
    def __init__(self, phi, e, w_s, v_sa, v_sd,  p_1, p_2) -> None:
        self.phi = phi
        self.e = e
        self.w_s = w_s
        self.v_sa = v_sa
        self.v_sd = v_sd

        self.p_1 =  p_1
        self.p_2 =  p_2


        '''
        Matriz de transformação da base A para a base B
            Base A: base ortonormal com origem em um vértice do quadrilátero.
                eixo x: apontondo para leste.
                eixo y: apontanto para o norte.
                eixo z: apontando para o zênite.

            Base B: rotação em (90 - phi)° da base A em torno do eixo x no sentido horário. Dessa forma o eixo z sempre está apontando para o polo Norte Celeste. Útil para descrever a posição do Sol.

        '''
        self.M_AB = M_AB = np.array([[1,0,0], [0,sin(phi), cos(phi)], [0 , -cos(phi), sin(phi)]])


    def d_t(self, t):
        '''
        Declinação em função do tempo (neagtivo quando o Sol está no hemisfério Sul Celeste).

        PARÂMETROS
        ----------
        t: float
            Instante do tempo.
        
        RETORNA
        -------
            Declinação do Sol no instante determinada.
        '''
        return arcsin(-cos(self.v_sa*t)*sin(self.e))


    # Ângulo horário (negativo a leste do meridiano local)
    def h_t(self, t):
        '''
        Ângulo horário em função do tempo (negativo a leste do meridiano local).

        PARÂMETROS
        ----------
        t: float
            Instante do tempo.
        
        RETORNA
        -------
            Ângulo horário do Sol no instante determinada.
        '''
        return -pi/2 + t*self.v_sd


    def area_t(self, t):
        '''
        Área do quadrilátero projetada em um plano ortogonal aos raios de luz em função do tempo
        
        PARÂMETROS
        ----------
        t: nparray
            Array contendo os instantes em que se deseja calcular a área.
            
        
        RETORNA
        -------
            A: 1-D nparray 
                Array com os valores das áreas para os instantes dados em t.
            s_t[2]: 1_D nparray
                Array com as posições z da base A (base no local do quadrilátero) para os intantes dados em t.
        '''

        s_tb = np.array([-cos(self.d_t(t))*sin(self.h_t(t)), -cos(self.d_t(t))*cos(self.h_t(t)) , sin(self.d_t(t))])
        s_t = self.M_AB.dot(s_tb)

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


        A = self.p_1.dot(v_1) * self.p_2.dot(v_2) - self.p_2.dot(v_1) * self.p_1.dot(v_2)
        return A, s_t[2]

    
    def time_interval_steps(self, delta_t, t_x_max, t_i, t_f):
        '''
        Seta o interalo de integração de area_t(t) e o tamnho de cada passo da integração

        PARÂMETROS
        ----------
            delta_t: float
                Comprimento dos retângulos utilizados na integração.
            t_x_max: int
                Número máximo de retângulos calculados a cada passo de integração.
            t_i: float
                Instante inicial.
            t_f: float
                Instante final.

        RETORNA
        -------
            Nada
        '''
        self.t_i = t_i # instante inicial
        self.t_f = t_f # instante final

        self.delta_t = delta_t # número máximo de retângulos calculados a cada passo de integração.
        self.size_step = (self.delta_t * t_x_max) # tamanho do intervalo de tempo de cada passo feito na integração.
        self.steps_t = np.floor( (t_f - t_i) / self.size_step ) # número de passos que serão necessários para realizar a integração (será necessário realizair um passo adicional se o intervalo de integração não for divisível por size_step)


    def calculate_energy(self):
        '''
        Calcula a enrgia solar que passa pelo quadrilátero no intervalo de tempo definido.

        PARÂMETROS
        ----------
            Nada
        
        RETORNA
        -------
            Nada
        '''
        int_area = 0
        t_x = [self.t_i]
        for _ in range(0, int(self.steps_t)):
            int_area_step = 0
            t_x = np.arange(t_x[-1], t_x[-1] + self.size_step, self.delta_t) + self.delta_t/2
            
            area_value, sun_z_pos = self.area_t(t_x)
            for area_i, s_z in zip(area_value, sun_z_pos):
                if s_z > 0:
                    if area_i > 0:
                        int_area_step += self.delta_t * area_i
            
            int_area += int_area_step

        if (self.size_step*self.steps_t < self.t_f):
            int_area_step = 0
            t_x = np.arange(t_x[-1], self.t_f, self.delta_t) + self.delta_t/2
            
            area_value, sun_z_pos = self.area_t(t_x)
            for area_i, s_z in zip(area_value, sun_z_pos):
                if s_z > 0:
                    if area_i > 0:
                        int_area_step += self.delta_t * area_i
            
            int_area += int_area_step

        total_energy = int_area * self.w_s * 24 * 60 * 60 * 1 / 3600000
        return total_energy


    def range_orientation(self, ang_1, ang_2, size_step):
        '''
        Calcular a energia incidente no quadrilátero para uma sequência de orientações. Quando ang_1 = ang_2 = 0, o quadrilátero está paralelo ao chão, com uma de suas arestas alinhada com o leste. Os resultados são printados no console.

        PARÂMETROS
        ----------
            and_1: float
                Ângulo entre a aresta que estava alinhada com o leste e a direção oeste-leste. Positivo no sentido anti-horário.
            and_2: float
                Inclinação do quadrilátero: ângulo da aresta que estava alinhada com o leste e o plano do horizonte. Positivo no sentido anti-horário.
            size_steps: int
                Diferença entre dois ângulos conssecutivos no intervalo de orientações.
        RETORNA
        -------
            Nada
        '''
        ang_1arr = np.arange(ang_1[0] , ang_1[1], size_step)
        ang_n1 = len(ang_1arr)
        ang_2arr = np.arange(ang_2[0], ang_2[1], size_step)
        ang_n2 = len(ang_2arr)

        ang_1arr, ang_2arr = np.meshgrid(ang_1arr, ang_2arr)


        energy = []
        for ang_1, ang_2 in zip(ang_1arr.reshape(ang_n1*ang_n2), ang_2arr.reshape(ang_n1*ang_n2)):
            ang_1, ang_2 = ang_1/180*pi, ang_2/180*pi

            M_pr = np.array([[cos(ang_1), -sin(ang_1), 0], [sin(ang_1), cos(ang_1), 0], [0,0,1]])
            self.p_1 = M_pr.dot(np.array([1,0,0]))
            self.p_2 = M_pr.dot(np.array([0, cos(ang_2), -sin(ang_2)]))

            total_energy = self.calculate_energy()

            print(f"{ang_1/pi*180};{ang_2/pi*180} | {total_energy}")
            energy.append(total_energy)

        ''' Gráfico dos resultados obtidos (futuramente criar um método dedica a essa função) '''
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

        energy = np.array(energy).reshape(ang_n2, ang_n1)
        
        # np.savetxt("ang_1.txt", ang_1arr)
        # np.savetxt("ang_2.txt", ang_2arr)
        # np.savetxt("energy.txt", energy)

        # print(ang_1arr.reshape(ang_n2, ang_n1))
        # print(ang_2arr.reshape(ang_n2, ang_n1))
        # print(energy)

        if len(energy) > 1:
            surf = ax.plot_surface(ang_1arr, ang_2arr, energy, cmap=cm.coolwarm)

            ax.set_xlabel("ang_1 (°)")
            ax.set_ylabel("ang_2 (°)")
            ax.set_zlabel("Energia (Kwh)")

            ax.zaxis.set_major_formatter('{x:.0f}')
            fig.colorbar(surf, shrink=0.5, aspect=5)

            plt.show()
