# simulation.py
import numpy as np

class PollutantModel:
    """
    Modelo de poluente em grade 2D usando esquema tipo autômato celular.
    Estados contínuos: concentração por célula (float >= 0).
    """
    def __init__(self, nx=200, ny=120, dx=1.0, dy=1.0, dt=1.0,
                 diffusion=0.1, decay=0.001, velocity_field=None,
                 boundary='open'):
        """
        nx, ny: dimensão da grade
        dx, dy: espaçamento espacial (usado para cálculo do laplaciano)
        dt: passo de tempo
        diffusion: coeficiente D (difusividade efetiva)
        decay: taxa lambda (por unidade de tempo)
        velocity_field: array shape (nx, ny, 2) ou function(i,j) -> (ux,uy)
        boundary: 'open', 'periodic', ou 'reflective'
        """
        self.nx = nx
        self.ny = ny
        self.dx = dx
        self.dy = dy
        self.dt = dt
        self.D = diffusion
        self.decay = decay
        self.boundary = boundary
        self.C = np.zeros((nx, ny), dtype=float)  # concentração
        if velocity_field is None:
            # default: uniform flow to the right
            self.u = np.zeros((nx, ny, 2), dtype=float)
            self.u[...,0] = 0.5  # ux
            self.u[...,1] = 0.0  # uy
        else:
            self.u = velocity_field if isinstance(velocity_field, np.ndarray) else self._build_velocity_from_func(velocity_field)
        # buffers para cálculos
        self._tmp = np.zeros_like(self.C)

    def _build_velocity_from_func(self, f):
        arr = np.zeros((self.nx, self.ny, 2), dtype=float)
        for i in range(self.nx):
            for j in range(self.ny):
                arr[i,j] = f(i,j)
        return arr

    def add_source_point(self, i, j, amount):
        """Adiciona instantaneamente amount (concentração) na célula (i,j)."""
        if 0 <= i < self.nx and 0 <= j < self.ny:
            self.C[i,j] += amount

    def _apply_boundary(self, arr):
        if self.boundary == 'periodic':
            # nada a fazer, índices modulo serão usados no cálculo
            return arr
        if self.boundary == 'reflective':
            # reflect edges by copying neighbor values
            arr[0,:] = arr[1,:]
            arr[-1,:] = arr[-2,:]
            arr[:,0] = arr[:,1]
            arr[:,-1] = arr[:,-2]
            return arr
        if self.boundary == 'open':
            # mantém zeros fora (flux sai da grade) -> já acontece se não "preenchemos"
            # para estabilidade podemos garantir borda = 0
            arr[0,:] = arr[0,:]
            arr[-1,:] = arr[-1,:]
            arr[:,0] = arr[:,0]
            arr[:,-1] = arr[:,-1]
            return arr
        return arr

    def step(self, source_func=None):
        """
        Executa um passo de tempo: advecção, difusão, decaimento e fontes.
        source_func: função opcional (t) -> lista de (i,j,amount) para adicionar
        """
        C = self.C
        nx, ny = self.nx, self.ny
        dt = self.dt
        dx, dy = self.dx, self.dy
        D = self.D
        u = self.u

        # ---------- Advecção (esquema flux-based simples) ----------
        # calculamos fluxos para 4 direções baseado em components u_x e u_y
        flux = np.zeros((nx, ny, 4), dtype=float)  # left, right, up, down
        # componente x
        ux = u[...,0]
        uy = u[...,1]

        # CFL-like fraction: fraction moved in dt = |u|*dt/dx (clamp to <=1)
        frac_x_pos = np.clip(ux * dt / dx, -1.0, 1.0)
        frac_y_pos = np.clip(uy * dt / dy, -1.0, 1.0)

        # fluxo para direita (from cell to right) quando ux>0
        # usamos aproximação: amount moved = C * max(0, frac_x_pos)
        move_right = C * np.maximum(frac_x_pos, 0.0)
        move_left  = C * np.maximum(-frac_x_pos, 0.0)
        move_down  = C * np.maximum(frac_y_pos, 0.0)
        move_up    = C * np.maximum(-frac_y_pos, 0.0)

        # inicia tmp com a parcela que permanece (resto após movimentos)
        remain = C - (move_right + move_left + move_down + move_up)

        # construir novo campo advectado
        C_adv = np.zeros_like(C)
        # cada célula recebe remanescentes
        C_adv += remain

        # deslocar move_right para a célula i+1
        C_adv[1:,:] += move_right[:-1,:]
        if self.boundary == 'periodic':
            C_adv[0,:] += move_right[-1,:]
        # move_left to i-1
        C_adv[:-1,:] += move_left[1:,:]
        if self.boundary == 'periodic':
            C_adv[-1,:] += move_left[0,:]
        # move_down to j+1
        C_adv[:,1:] += move_down[:,:-1]
        if self.boundary == 'periodic':
            C_adv[:,0] += move_down[:,-1]
        # move_up to j-1
        C_adv[:,:-1] += move_up[:,1:]
        if self.boundary == 'periodic':
            C_adv[:,-1] += move_up[:,0]

        # para fronteiras 'open' a massa movida para fora é perdida (simula escoamento para longe)
        # para 'reflective' poderíamos devolver massa — omitido para simplicidade

        # ---------- Difusão (laplaciano discreto) ----------
        lap = np.zeros_like(C_adv)
        # vizinhança 4 (von Neumann)
        lap[1:-1,1:-1] = (C_adv[2:,1:-1] + C_adv[:-2,1:-1] +
                         C_adv[1:-1,2:] + C_adv[1:-1,:-2] -
                         4*C_adv[1:-1,1:-1]) / (dx*dy)
        # bordas (trata por condições de contorno simples)
        if self.boundary == 'periodic':
            # computa laplaciano com periodicidade (simplificado)
            lap = (np.roll(C_adv,1,axis=0) + np.roll(C_adv,-1,axis=0) +
                   np.roll(C_adv,1,axis=1) + np.roll(C_adv,-1,axis=1) -
                   4*C_adv) / (dx*dy)
        # termo difusivo
        C_diff = C_adv + D * dt * lap

        # ---------- Decaimento ----------
        C_decay = C_diff * np.exp(-self.decay * dt)

        # ---------- Fonte (opcional) ----------
        if source_func is not None:
            items = source_func()
            for (i,j,amt) in items:
                if 0 <= i < nx and 0 <= j < ny:
                    C_decay[i,j] += amt

        # garantir não-negatividade
        C_decay[C_decay < 0] = 0.0

        # aplicar condição de contorno
        C_final = self._apply_boundary(C_decay)

        self.C = C_final
        return self.C

    def total_mass(self):
        return np.sum(self.C) * self.dx * self.dy
