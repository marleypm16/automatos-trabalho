# Simulação 2D de Propagação de Poluentes

Este projeto implementa um modelo de transporte e difusão de poluentes em uma grade 2D, com visualização interativa em tempo real e geração de gráficos para análise dos resultados.

## 📋 Requisitos

- Python 3.9 ou superior

### Dependências

Instale os pacotes necessários (recomenda-se usar um ambiente virtual):
```bash
pip install numpy matplotlib pygame
```

## 🗂️ Estrutura do Projeto
```
.
├── simulation.py           # Classe PollutantModel (lógica da simulação)
├── visualize_pygame.py     # Interface interativa com pygame
└── plots.py                # Geração de gráficos estáticos
```

### Descrição dos Módulos

- **`simulation.py`**: Implementa a classe `PollutantModel` com os processos físicos de advecção, difusão, decaimento e emissão de fontes poluidoras.

- **`visualize_pygame.py`**: Interface gráfica interativa que executa a simulação em tempo real, exibe o campo de concentração e permite salvar o estado final.

- **`plots.py`**: Gera visualizações estáticas a partir dos dados salvos:
  - `history_mass.npy`: evolução da massa total ao longo do tempo
  - `final_concentration.npy`: distribuição espacial final da concentração

## 🚀 Como Usar

### 1. Executar a simulação interativa
```bash
python visualize_pygame.py
```

A simulação será exibida em uma janela do pygame mostrando a propagação dos poluentes em tempo real.

### 2. Gerar os gráficos de análise

Após executar a simulação interativa, dois arquivos serão criados automaticamente (`history_mass.npy` e `final_concentration.npy`). Para visualizar os gráficos:
```bash
python plots.py
```

Isso gerará visualizações estáticas da evolução temporal da massa total e mapa de calor da concentração final.
