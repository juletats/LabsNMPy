def findIndex(x, massive):
    index = 0
    for i in range(len(massive)):
        if x == massive[i]:
            index = i;
    return index


x = [-1.0, 0.0, 1.0, 2.0, 3.0]
y = [2.3562, 1.5708, 0.7854, 0.46365, 0.32175]
X = 1

index = findIndex(X, x)

leftDer_y = (y[index] - y[index - 1]) / (x[index] - x[index - 1])  # левая производная
rightDer_y = (y[index + 1] - y[index]) / (x[index + 1] - x[index])  # правая производная

Der_y = leftDer_y + (((rightDer_y - leftDer_y) / (x[index + 1] - x[index - 1]))) * (
            2 * X - x[index - 1] - x[index])  # сама проивзодная
print("\nпервая производная:       y' = ", Der_y)

DerSecond_y = 2 * ((rightDer_y - leftDer_y) / (x[index + 1] - x[index - 1]))  # вторая проивзодная
print("\nпервая производная:      y'' = ", DerSecond_y)

print()
