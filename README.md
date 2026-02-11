# Autômato Celular Epidemiológico

Simulação de propagação epidêmica usando autômatos celulares com modelo SID (Suscetível-Infectado-Morto).

## Descrição

Este projeto implementa um autômato celular 2D que simula a dinâmica de uma epidemia em uma população. O modelo considera três estados:

- **S (Suscetível)**: Indivíduos que podem ser infectados (azul)
- **I (Infectado)**: Indivíduos infectados que podem infectar vizinhos (vermelho)
- **D (Morto)**: Indivíduos que morreram pela doença (preto)

## Funcionalidades

O código gera:

1. **Curvas de Entropia H(t)**: Mede a desordem do sistema ao longo do tempo para diferentes taxas de mortalidade (ρ)
2. **Curvas de Prevalência**: Evolução temporal da porcentagem de infectados
3. **Curvas de Óbitos**: Evolução temporal da mortalidade acumulada
4. **Snapshots Espaciais**: Visualização da grade em diferentes momentos (t=50, 100, 200) para três regimes de mortalidade
5. **Arquivo de Resultados**: Dados numéricos salvos em `resultados_simulacao.txt`

## Parâmetros da Simulação

- **L**: Tamanho da grade (100×100)
- **T**: Número de passos de tempo (300)
- **β (beta)**: Taxa de transmissão (0.3)
- **ρ (rho)**: Taxa de mortalidade - testado com três valores:
  - 0.0 (baixa mortalidade)
  - 0.3 (mortalidade média)
  - 0.7 (alta mortalidade)

## Como Rodar

### Pré-requisitos

- Python 3.7 ou superior
- pip (gerenciador de pacotes Python)

### Instalação das Dependências

```bash
pip install -r requirements.txt
```

Ou instale manualmente:

```bash
pip install numpy matplotlib
```

### Execução

```bash
python main.py
```

## Saída

Após executar o código, você obterá:

1. **3 figuras interativas** mostrando:
   - Entropia de Shannon ao longo do tempo
   - Prevalência de infectados e óbitos acumulados
   - Padrões espaciais da grade em diferentes momentos

2. **Arquivo `resultados_simulacao.txt`** contendo:
   - Estatísticas resumidas para cada valor de ρ
   - Estado final da população
   - Pico de infectados
   - Séries temporais com amostragem a cada 10 passos

## Modelo Matemático

### Regras de Transição

- **S → I**: Suscetível infecta com probabilidade $P = 1 - (1-\beta)^k$, onde k é o número de vizinhos infectados (vizinhança de Moore)
- **I → D**: Infectado morre com probabilidade ρ a cada passo de tempo
- **D → D**: Mortos permanecem mortos

### Contorno Periódico

A grade usa contorno periódico (toroidal), onde as bordas se conectam.

## Resultados Esperados

- **ρ = 0.0**: Epidemia se espalha por toda população sem óbitos
- **ρ = 0.3**: Epidemia rápida com alta mortalidade (~99%)
- **ρ = 0.7**: Mortalidade muito alta limita a propagação, mas ainda atinge ~97% de óbitos
