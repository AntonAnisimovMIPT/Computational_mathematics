import matplotlib.pyplot as plt

x_values = []
y_values = []

with open("data.txt", "r") as file:
    for line in file:
        x, y = map(float, line.split())
        x_values.append(x)
        y_values.append(y)

plt.plot(x_values, y_values, label="Решение", color="b")

plt.xlabel("x")
plt.ylabel("y")
plt.title("График решения")
plt.legend()

plt.savefig("plot.png", dpi=300)

plt.close()
