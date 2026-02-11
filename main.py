import numpy as np
import matplotlib.pyplot as plt

# Estados: 0 = morto (D), 1 = suscetível (S), 2 = infectado (I)


def inicializar_grade(L, p_inicial_infectados=0.01, seed=None):
    if seed is not None:
        np.random.seed(seed)
    grade = np.ones((L, L), dtype=np.int8)  # todo mundo suscetível
    mascara_infectados = np.random.rand(L, L) < p_inicial_infectados
    grade[mascara_infectados] = 2
    return grade


def vizinhos_infectados(grade):
    """Conta vizinhos infectados (Moore, contorno periódico)."""
    I = (grade == 2).astype(np.int8)
    # shifts em 8 direções
    neigh_sum = (
        np.roll(I, 1, axis=0)
        + np.roll(I, -1, axis=0)
        + np.roll(I, 1, axis=1)
        + np.roll(I, -1, axis=1)
        + np.roll(np.roll(I, 1, axis=0), 1, axis=1)
        + np.roll(np.roll(I, 1, axis=0), -1, axis=1)
        + np.roll(np.roll(I, -1, axis=0), 1, axis=1)
        + np.roll(np.roll(I, -1, axis=0), -1, axis=1)
    )
    return neigh_sum


def passo(grade, beta, rho):
    L = grade.shape[0]
    nova = grade.copy()

    # Máscaras de estado atual
    S = grade == 1
    I = grade == 2
    D = grade == 0

    # 1) Suscetíveis podem ser infectados
    k = vizinhos_infectados(grade)
    # Probabilidade de não infecção por cada vizinho: (1 - beta)^k
    p_infeccao = 1.0 - (1.0 - beta) ** k
    rand_inf = np.random.rand(L, L)
    infectam = S & (rand_inf < p_infeccao)

    # 2) Infectados podem morrer
    rand_morte = np.random.rand(L, L)
    morrem = I & (rand_morte < rho)

    # Atualiza estados
    nova[infectam] = 2
    nova[morrem] = 0
    # Infectados que não morreram permanecem infectados (já estão em 2)
    # Mortos permanecem mortos; suscetíveis que não infectam permanecem em 1

    return nova


def entropia_shannon(grade):
    L2 = grade.size
    pS = np.sum(grade == 1) / L2
    pI = np.sum(grade == 2) / L2
    pD = np.sum(grade == 0) / L2
    probs = np.array([pS, pI, pD])
    probs = probs[probs > 0.0]
    H = -np.sum(probs * np.log2(probs))
    return H, (pS, pI, pD)


def simular(
    L=100,
    T=300,
    beta=0.3,
    rho=0.3,
    p_inicial_infectados=0.01,
    seed=None,
    snapshots_tempos=None,
):
    grade = inicializar_grade(L, p_inicial_infectados, seed)
    Hs = []
    fracS = []
    fracI = []
    fracD = []
    snapshots = {}

    for t in range(T):
        H, (pS, pI, pD) = entropia_shannon(grade)
        Hs.append(H)
        fracS.append(pS)
        fracI.append(pI)
        fracD.append(pD)

        # Capturar snapshots em tempos específicos
        if snapshots_tempos and t in snapshots_tempos:
            snapshots[t] = grade.copy()

        grade = passo(grade, beta, rho)

    return np.array(Hs), np.array(fracS), np.array(fracI), np.array(fracD), snapshots


def salvar_resultados(filename, rhos, resultados):
    """Salva os resultados das simulações em arquivo .txt"""
    with open(filename, "w", encoding="utf-8") as f:
        f.write("=" * 80 + "\n")
        f.write("RESULTADOS DA SIMULAÇÃO - AUTÔMATO CELULAR EPIDEMIOLÓGICO\n")
        f.write("=" * 80 + "\n\n")

        for i, rho in enumerate(rhos):
            Hs, fS, fI, fD, snapshots = resultados[i]

            f.write(f"\n{'='*80}\n")
            f.write(f"PARÂMETRO: ρ (mortalidade) = {rho}\n")
            f.write(f"{'='*80}\n\n")

            # Estatísticas resumidas
            f.write(f"--- Estatísticas da Entropia H(t) ---\n")
            f.write(f"  H inicial:  {Hs[0]:.6f}\n")
            f.write(f"  H final:    {Hs[-1]:.6f}\n")
            f.write(f"  H máxima:   {np.max(Hs):.6f}\n")
            f.write(f"  H mínima:   {np.min(Hs):.6f}\n")
            f.write(f"  H média:    {np.mean(Hs):.6f}\n\n")

            f.write(f"--- Estado Final da População (t={len(Hs)-1}) ---\n")
            f.write(f"  Suscetíveis (S): {fS[-1]*100:.2f}%\n")
            f.write(f"  Infectados (I):  {fI[-1]*100:.2f}%\n")
            f.write(f"  Mortos (D):      {fD[-1]*100:.2f}%\n\n")

            f.write(f"--- Pico de Infectados ---\n")
            pico_I_idx = np.argmax(fI)
            f.write(f"  Tempo do pico:     t = {pico_I_idx}\n")
            f.write(f"  Prevalência máx:   {fI[pico_I_idx]*100:.2f}%\n\n")

            f.write(f"--- Total de Óbitos Acumulados ---\n")
            f.write(f"  Mortalidade total: {fD[-1]*100:.2f}%\n\n")

            # Séries temporais (amostragem a cada 10 passos para não ficar muito grande)
            f.write(f"--- Séries Temporais (amostragem a cada 10 passos) ---\n")
            f.write(
                f"{'Tempo':>6} | {'H(t)':>10} | {'S(%)':>8} | {'I(%)':>8} | {'D(%)':>8}\n"
            )
            f.write(f"{'-'*6}-+-{'-'*10}-+-{'-'*8}-+-{'-'*8}-+-{'-'*8}\n")
            for t in range(0, len(Hs), 10):
                f.write(
                    f"{t:6d} | {Hs[t]:10.6f} | {fS[t]*100:8.2f} | {fI[t]*100:8.2f} | {fD[t]*100:8.2f}\n"
                )

            f.write(f"\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("FIM DO RELATÓRIO\n")
        f.write("=" * 80 + "\n")


if __name__ == "__main__":
    L = 100
    T = 300
    beta = 0.3

    # Comparar três valores de mortalidade
    rhos = [0.0, 0.3, 0.7]
    cores = ["tab:blue", "tab:orange", "tab:green"]
    labels = ["Baixa (ρ=0.0)", "Média (ρ=0.3)", "Alta (ρ=0.7)"]

    # Tempos para capturar snapshots da grade
    snapshots_tempos = [50, 100, 200]

    # Executar simulações
    print("Executando simulações...")
    resultados = []
    for rho in rhos:
        print(f"  Simulando com ρ = {rho}...")
        Hs, fS, fI, fD, snapshots = simular(
            L=L,
            T=T,
            beta=beta,
            rho=rho,
            p_inicial_infectados=0.01,
            seed=42,
            snapshots_tempos=snapshots_tempos,
        )
        resultados.append((Hs, fS, fI, fD, snapshots))

    # Salvar resultados em arquivo .txt
    print("\nSalvando resultados em 'resultados_simulacao.txt'...")
    salvar_resultados("resultados_simulacao.txt", rhos, resultados)
    print("Arquivo salvo com sucesso!")

    # FIGURA 1: Entropia H(t)
    plt.figure(figsize=(10, 4))
    for i, (rho, cor, label) in enumerate(zip(rhos, cores, labels)):
        Hs = resultados[i][0]
        plt.plot(Hs, label=label, color=cor, linewidth=2)
    plt.xlabel("Tempo (passos)", fontsize=12)
    plt.ylabel("Entropia de Shannon H(t)", fontsize=12)
    plt.title("Evolução da Entropia para Diferentes Mortalidades", fontsize=13)
    plt.legend(fontsize=10)
    plt.grid(alpha=0.3)
    plt.tight_layout()

    # FIGURA 2: Prevalência (Infectados) e Óbitos
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Subplot 1: Prevalência de infectados
    for i, (rho, cor, label) in enumerate(zip(rhos, cores, labels)):
        fI = resultados[i][2]
        ax1.plot(fI * 100, label=label, color=cor, linewidth=2)
    ax1.set_xlabel("Tempo (passos)", fontsize=12)
    ax1.set_ylabel("Prevalência de Infectados (%)", fontsize=12)
    ax1.set_title("Curvas de Prevalência (I)", fontsize=13)
    ax1.legend(fontsize=10)
    ax1.grid(alpha=0.3)

    # Subplot 2: Óbitos acumulados
    for i, (rho, cor, label) in enumerate(zip(rhos, cores, labels)):
        fD = resultados[i][3]
        ax2.plot(fD * 100, label=label, color=cor, linewidth=2)
    ax2.set_xlabel("Tempo (passos)", fontsize=12)
    ax2.set_ylabel("Óbitos Acumulados (%)", fontsize=12)
    ax2.set_title("Curvas de Óbitos (D)", fontsize=13)
    ax2.legend(fontsize=10)
    ax2.grid(alpha=0.3)

    plt.tight_layout()

    # FIGURA 3: Snapshots espaciais da grade
    fig, axes = plt.subplots(len(rhos), len(snapshots_tempos), figsize=(12, 10))

    cmap = plt.cm.colors.ListedColormap(["black", "lightblue", "red"])
    bounds = [-0.5, 0.5, 1.5, 2.5]
    norm = plt.cm.colors.BoundaryNorm(bounds, cmap.N)

    for i, (rho, label) in enumerate(zip(rhos, labels)):
        snapshots = resultados[i][4]
        for j, t in enumerate(snapshots_tempos):
            ax = axes[i, j] if len(rhos) > 1 else axes[j]
            im = ax.imshow(snapshots[t], cmap=cmap, norm=norm, interpolation="nearest")

            # Título das colunas
            if i == 0:
                ax.set_title(f"t = {t}", fontsize=11)

            # Rótulo das linhas
            if j == 0:
                ax.set_ylabel(label, fontsize=11, rotation=90, labelpad=10)

            ax.set_xticks([])
            ax.set_yticks([])

    # Adicionar legenda de cores
    fig.suptitle(
        "Padrões Espaciais da Grade (Preto=Morto, Azul=Suscetível, Vermelho=Infectado)",
        fontsize=13,
        y=0.98,
    )
    plt.tight_layout()

    print("\nExibindo gráficos...")
    plt.show()

    print("\n✓ Simulação concluída!")
    print(f"✓ Resultados salvos em: resultados_simulacao.txt")
