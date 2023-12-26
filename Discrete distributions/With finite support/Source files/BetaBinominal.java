public class BetaBinominal extends MathLibrary {
    int n;
    double alpha;
    double beta;

    BetaBinominal(int n, double alpha, double beta) {
        if (n <= 0) throw new IllegalArgumentException("n > 0");
        if (alpha < 0) throw new IllegalArgumentException("alpha > 0");
        if (beta < 0) throw new IllegalArgumentException("beta > 0");
        this.n = n;
        this.alpha = alpha;
        this.beta = beta;
    }

    double PMF(int x) {
        if (x < 0 || x > n) throw new IllegalArgumentException("0 <= x <= n");
        return binomial(n, x) * (beta(x + alpha, n - x + alpha) / beta(alpha, beta));
    }

    double CDF(int x) {
        if (x < 0 || x > n) throw new IllegalArgumentException("0 <= x <= n");
        double value = 0;
        for (int i = 0; i < x; i++) {
            value += (gamma(alpha + i) * gamma(beta + n - i)) / (gamma(i + 1) * gamma(n - i + 1));
        }
        return (1 / beta(alpha, beta)) * (gamma(n + 1) / (gamma(alpha + beta + n))) * value;
    }

    double Mean() {
        return (n * alpha) / (alpha + beta);
    }

    double Variance() {
        double licznik = n * alpha * beta * (alpha + beta + n);
        double mianownik = power(alpha + beta, 2) * (alpha + beta + 1);
        return licznik / mianownik;
    }

    double Skewness() {
        double licznik1 = (alpha + beta + (2 * n)) * (beta - alpha);
        double wartosc1 = licznik1 / (alpha + beta + 2);
        double mianownik2 = n * alpha * beta * (n + alpha + beta);
        double wartosc2 = root((1 * alpha * beta) / mianownik2, 2);
        return wartosc1 * wartosc2;
    }

    double MGF(double t) {
        return hypergeometric(-1 * n, alpha, alpha + beta, 1 - euler(t));
    }

    double PGF(double z) {
        return hypergeometric(-1 * n, alpha, alpha + beta, 1 - z);
    }
}
