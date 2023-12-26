public class FishersNoncentralHypergeometric extends MathLibrary{
    int m1;
    int m2;
    int N;
    int n;
    double omega;
    int min;
    int max;

    FishersNoncentralHypergeometric(int m1, int m2, int n, double omega) {
        if (m1 < 0 || m2 < 0) throw new IllegalArgumentException("m1, m2 >= 0");
        if (n < 0 || n > (m1 + m2)) throw new IllegalArgumentException("0 <= n <= m1 + m2");
        if (omega < 0) throw new IllegalArgumentException("omega >= 0");
        this.m1 = m1;
        this.m2 = m2;
        this.N = m1 + m2;
        this.n = n;
        this.omega = omega;
        this.min = (int) max(0, n - m2);
        this.max = (int) min(n, m1);
    }

    private double support(int k) {
        double wartosc = 0;
        for (int i = min; i < max; i++) {
            wartosc += binomial(m1, i) * binomial(m2, n - i) * power(omega, i) * power(i, k);
        }
        return wartosc;
    }

    double PMF(int x) {
        return (binomial(m1, x) * binomial(m2, n - x) * power(omega, x)) / support(0);
    }

    double Mean() {
        return support(1) / support(0);
    }

    int Mode() {
        double A = omega - 1;
        double B = m1 + n - N -1 * (m1 + n + 2) * omega;
        double C = (m1 + 1) * (n + 1) * omega;
        return (int) ((-2 * C) / (B - root(power(B, 2) - 4 * A * C, 2)));
    }

    double Variance() {
        return (support(2) / support(0)) - power(support(1) / support(0), 2);
    }
}
