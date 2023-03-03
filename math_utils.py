import numpy as np
from numpy import cos, sin

def rotation_matrix(rotation: list[tuple]) -> np.ndarray:
    '''
    Retorna a matrix que contém as rotações em série em "rotation".

    Parâmetros:
        rotation: list[tuples]
            Cada elemento dessa lista é uma tupla cujo primeiro elemento é o eixo de rotação
            e o segundo elemento o ângulo de rotação em radianos.
            O i-ésimo elemento dessa lista contém a i-ésima rotação que é aplicada.
    '''
    rot_final = np.identity(3)
    for axis, angle in rotation:
        if axis == "x":
            rot_el = np.array([[1, 0, 0], [0, cos(angle), -sin(angle)], [0, sin(angle), cos(angle)]])
        elif axis == "y":
            rot_el = np.array([[cos(angle), 0, sin(angle)], [0, 1, 0], [-sin(angle), 0, cos(angle)]])
        elif axis == "z":
            rot_el = np.array([[cos(angle), -sin(angle), 0], [sin(angle), cos(angle), 0], [0,0, 1]])

        rot_final = rot_el.dot(rot_final)

    return rot_final

def numeric_integral(f, dx):
    '''  
    Sendo f uma lista de pontos de uma função, cujo o espaçamento 
    entre dois pontos é "dx", é calculado a integral dessa função, 
    utilizando o método dos trapézios, para o intervalo dos pontos dados.
    '''
    return ((f[0] + f[-1])/2 + sum(f[1:-1])) * dx

