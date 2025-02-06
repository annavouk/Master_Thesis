import prince
import pandas as pd
import matplotlib.pyplot as plt
import seaborn

df = pd.read_pickle('seeds_binary_per_patric.pckl')

mca = prince.MCA(n_components=2, random_state=42)
mca = mca.fit(df)

mca_results = mca.transform(df)

eigenvalues = mca.eigenvalues_
explained_inertia = eigenvalues / sum(eigenvalues)

print(mca_results.head())

plt.figure(figsize=(12, 8))

plt.scatter(
    mca_results.iloc[:, 0], 
    mca_results.iloc[:, 1], 
    c="blue", alpha=0.6, edgecolor="k", s=50
)

plt.title('MCA - Multiple Correspondence Analysis', fontsize=16)
plt.xlabel(f'Component 1 ({explained_inertia[0]*100:.2f}% variance)', fontsize=14)
plt.ylabel(f'Component 2 ({explained_inertia[1]*100:.2f}% variance)', fontsize=14)

plt.axhline(0, color="black", linewidth=0.5, linestyle="--")
plt.axvline(0, color="black", linewidth=0.5, linestyle="--")
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()

correlation_matrix = mca_results.corr()

plt.figure(figsize=(8, 6))
sns.heatmap(correlation_matrix, annot=True, cmap="coolwarm", fmt=".2f", linewidths=0.5)

plt.title('Heatmap of Component Correlations', fontsize=16)
plt.tight_layout()
plt.show()
