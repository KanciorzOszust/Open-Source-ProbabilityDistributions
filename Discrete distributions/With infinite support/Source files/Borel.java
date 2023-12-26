public class Borel extends MathLibrary {
    double mu;
    Borel(double mu) {
        if (mu < 0 || mu > 1) throw new IllegalArgumentException("0 <= m <= 1");
        this.mu = mu;
    }

    double PMF(int n) {
        double licznik = euler(-1 * mu * n) * power(mu * n, n - 1);
        return licznik / factorial(n);
    }

    double Mean() {
        return 1 / (1 - mu);
    }

    double Variance() {
        return mu / power(1 - mu, 3);
    }
}
