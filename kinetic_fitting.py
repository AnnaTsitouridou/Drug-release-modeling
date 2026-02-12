import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

# Make sure x and y are numpy arrays and flattened
x = np.array(x).flatten()
y = np.array(y).flatten()

xfit = np.linspace(min(x), max(x), 300)

# --------------------------------------------------
# 1. Define kinetic model equations
# --------------------------------------------------

# Zero-order
def zero_order(x, k0):
    return k0 * x

# First-order
def first_order(x, Cinf, k):
    return Cinf * (1 - np.exp(-k * x))

# Higuchi
def higuchi(x, kH):
    return kH * np.sqrt(x)

# Korsmeyer-Peppas
def korsmeyer_peppas(x, kKP, n):
    return kKP * (x ** n)

# --------------------------------------------------
# 2. Fit each model
# --------------------------------------------------

popt_zero, _ = curve_fit(zero_order, x, y)
popt_first, _ = curve_fit(first_order, x, y, p0=[max(y), 0.1])
popt_higuchi, _ = curve_fit(higuchi, x, y)
popt_kp, _ = curve_fit(korsmeyer_peppas, x, y)

# Predictions
y0 = zero_order(x, *popt_zero)
y1 = first_order(x, *popt_first)
yh = higuchi(x, *popt_higuchi)
ykp = korsmeyer_peppas(x, *popt_kp)

# --------------------------------------------------
# 3. Plot all fits
# --------------------------------------------------

plt.figure()
plt.scatter(x, y, label="Data")

plt.plot(xfit, zero_order(xfit, *popt_zero), label="Zero-order")
plt.plot(xfit, first_order(xfit, *popt_first), label="First-order")
plt.plot(xfit, higuchi(xfit, *popt_higuchi), label="Higuchi")
plt.plot(xfit, korsmeyer_peppas(xfit, *popt_kp), label="Korsmeyer-Peppas")

plt.xlabel("Time")
plt.ylabel("Fraction Released")
plt.title("Drug Release Kinetic Model Fitting")
plt.legend()
plt.grid()
plt.show()

# --------------------------------------------------
# 4. Calculate R² and RMSE
# --------------------------------------------------

def R_squared(y_true, y_pred):
    return 1 - np.sum((y_true - y_pred)**2) / np.sum((y_true - np.mean(y_true))**2)

rmse0  = np.sqrt(np.mean((y - y0)**2))
rmse1  = np.sqrt(np.mean((y - y1)**2))
rmseh  = np.sqrt(np.mean((y - yh)**2))
rmseKP = np.sqrt(np.mean((y - ykp)**2))

R2_0  = R_squared(y, y0)
R2_1  = R_squared(y, y1)
R2_H  = R_squared(y, yh)
R2_KP = R_squared(y, ykp)

# --------------------------------------------------
# 5. Build comparison table
# --------------------------------------------------

ModelNames = ["Zero-order", "First-order", "Higuchi", "Korsmeyer-Peppas"]

Params = [
    f"k0={popt_zero[0]:.4f}",
    f"Cinf={popt_first[0]:.4f}, k={popt_first[1]:.4f}",
    f"kH={popt_higuchi[0]:.4f}",
    f"kKP={popt_kp[0]:.4f}, n={popt_kp[1]:.4f}"
]

R2_values = [R2_0, R2_1, R2_H, R2_KP]
RMSE_values = [rmse0, rmse1, rmseh, rmseKP]

ResultsTable = pd.DataFrame({
    "Model": ModelNames,
    "Parameters": Params,
    "R²": R2_values,
    "RMSE": RMSE_values
})

print("-------------------------------------------------------")
print(ResultsTable)
print("-------------------------------------------------------")
print("Lower RMSE and higher R² indicate a better fit.")
