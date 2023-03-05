# painel-solar
Script para calcular quanta energia solar é incidida em uma planta solar, em qualquer lugar do planeta e em qualquer período de tempo.

Sete as constantes, a orientação e quantidade das placas solares e o intervalo de tempo no arquivo main.py e rode o mesmo para obter o resultado.

O ponto inicial do tempo é quando o Sol está com máxima declinação (solstício de inverno no hemisfério Sul e solstício de verão no hemisfério Norte).

## Pressuposições: 
* É assumido clima perfeito (a luz do Sol nunca é bloqueada por fatores climáticos).
* Nenhum objeto está entre as placas solares e o Sol durante todo o ano.
* Considerado velocidade angular constante para a translação e rotação da Terra.
* Não é considerado a variação da irradiação solar devido a excentricidade da órbita da Terra.
* É utilizado a seguinte expressão para o cálculo da irradiação solar.
$$I = 1,353~~kW/m^2 ~\cdot 0,7^{AM^{0,678}}$$
Em que $AM$ é a massa de ar, dada por
$$AM = \frac{1}{cos(z_s)}$$
Em que $z_s$ é o ângulo zenital do Sol.


