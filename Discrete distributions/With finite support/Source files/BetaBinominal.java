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
        return binomial(n, x) * (beta(x + alpha, n - x + alpha) / beta(alpha, beta));
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
}
