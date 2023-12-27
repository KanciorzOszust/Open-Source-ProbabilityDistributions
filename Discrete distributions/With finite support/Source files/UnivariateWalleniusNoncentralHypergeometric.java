public class UnivariateWalleniusNoncentralHypergeometric extends MathLibrary{
    int m1;
    int m2;
    int N;
    double n;
    double omega;

    UnivariateWalleniusNoncentralHypergeometric(int m1, int m2, double n, int omega) {
        if (n < 0 || n >= (m1 + m2)) throw new IllegalArgumentException("0 <= n < m1 + m2");
        if (omega <= 0) throw new IllegalArgumentException("omega > 0");
        this.m1 = m1;
        this.m2 = m2;
        this.N = m1 + m2;
        this.n = n;
        this.omega = omega;
    }

    private double D(double x) {
        return omega * (m1 - x) + (m2 - (n - x));
    }

    private double integral(double x) {
        double stepsize = 1.0 / 10000;
        double value = 0;
        for (double t = 0; t <= 1; t += stepsize) {
            value += doublePower(1 - doublePower(t, omega / D(x)), x) * doublePower(1 - doublePower(t, 1 / D(x)), n - x);
        }
        return value;
    }

    double PMF(double x) {
        if (x < max(0, n - m2) || x > min(n, m1)) {
            throw new IllegalArgumentException("max(0, n - m2) <= x <= min(n, m1)");
        }
        return binomial(m1, x) * binomial(m2, n - x) * integral(x);
    }
}
