# visualize_pygame.py
import sys
import pygame
import numpy as np
from simulation import PollutantModel
import matplotlib.cm as cm

def concentration_to_surface(C, width=800, height=480, cmap_name='viridis'):
    """
    Converte matriz de concentração C (nx,ny) para superfície pygame (scale to size).
    """
    nx, ny = C.shape
    # normalizar entre 0 e 1 para colormap
    vmin = 0.0
    vmax = max(1e-6, C.max())
    norm = (C - vmin) / (vmax - vmin)
    cmap = cm.get_cmap(cmap_name)
    rgba = (cmap(norm) * 255).astype(np.uint8)  # shape (nx,ny,4)
    # criar superfície
    surf = pygame.surfarray.make_surface(np.flipud(np.transpose(rgba[...,:3], (1,0,2))))
    surf = pygame.transform.scale(surf, (width, height))
    return surf

def main():
    pygame.init()
    width, height = 900, 540
    screen = pygame.display.set_mode((width, height))
    pygame.display.set_caption("Simulação: Propagação de Poluentes (AC)")

    # parâmetros do modelo
    nx, ny = 180, 110
    model = PollutantModel(nx=nx, ny=ny, dt=1.0, diffusion=0.2, decay=0.0010, boundary='open')

    # Exemplo de campo de velocidade: uniforme para direita com um pequeno rotacional no centro
    u = np.zeros((nx, ny, 2), dtype=float)
    u[...,0] = 0.8  # ux
    u[...,1] = 0.0
    # cria um "eddy" (vórtice) no centro
    cx, cy = nx//2, ny//2
    for i in range(nx):
        for j in range(ny):
            dx = i - cx
            dy = j - cy
            r2 = dx*dx + dy*dy + 1e-6
            # rotacional local decrescente com r
            u[i,j,0] += -0.6 * dy / np.sqrt(r2) * np.exp(-r2/(2*(min(nx,ny)/6)**2))
            u[i,j,1] +=  0.6 * dx / np.sqrt(r2) * np.exp(-r2/(2*(min(nx,ny)/6)**2))
    model.u = u

    clock = pygame.time.Clock()
    paused = False
    t = 0
    history = []  # salva massa total ao longo do tempo

    # adiciona fonte pontual constante na borda esquerda (simula descarga)
    def source():
        # devolve uma lista de (i,j,amount)
        out = []
        for jj in range(ny//2 - 2, ny//2 + 3):
            out.append((5, jj, 1.5))  # ajustável
        return out

    # inicial: pulso inicial no canto esquerdo
    for j in range(ny//2-1, ny//2+2):
        model.add_source_point(3, j, 50.0)

    running = True
    font = pygame.font.SysFont(None, 22)
    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_SPACE:
                    paused = not paused
                elif event.key == pygame.K_q:
                    running = False
                elif event.key == pygame.K_s:
                    # salva frame atual como imagem
                    pygame.image.save(screen, f"frame_{t:05d}.png")
        if not paused:
            C = model.step(source_func=source)
            history.append(model.total_mass())
            t += 1

        surf = concentration_to_surface(model.C, width=width, height=height, cmap_name='viridis')
        screen.blit(surf, (0,0))

        info = f"t={t}  total_mass={history[-1]:.2f}  max={model.C.max():.2f}  paused={'ON' if paused else 'OFF'}"
        txt = font.render(info, True, (255,255,255))
        screen.blit(txt, (5,5))

        pygame.display.flip()
        clock.tick(30)  # limitar a 30 FPS

    pygame.quit()
    # ao finalizar, salva histórico
    np.save("history_mass.npy", np.array(history))
    np.save("final_concentration.npy", model.C)
    print("Simulação encerrada. Arquivos salvos: history_mass.npy, final_concentration.npy")

if __name__ == "__main__":
    main()
