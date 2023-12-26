public class BetaNegativeBinomial extends MathLibrary {
    double alpha;
    double beta;
    int r;

    BetaNegativeBinomial(double alpha, double beta, int r) {
        if (alpha <= 0) throw new IllegalArgumentException("alpha > 0");
        if (beta <= 0) throw new IllegalArgumentException("beta > 0");
        if (r <= 0) throw new IllegalArgumentException("r > 0");
        this.alpha = alpha;
        this.beta = beta;
        this.r = r;
    }

    double PMF(int k) {
        double part1 = beta(r + k, alpha + beta) / (double) beta(r, alpha);
        double part2 = (double) gamma(k + beta) / (factorial(k) * gamma(beta));
        return part1 * part2;
    }

    double Mean() {
        if (alpha > 1) {
            return (double) (r * beta) / (alpha - 1);
        } else {
            return Double.POSITIVE_INFINITY;
        }
    }

    double Variance() {
        if (alpha > 2) {
            double licznik = r * beta * (r + alpha - 1) * (beta + alpha - 1);
            double mianownik = (alpha - 2) * power(alpha - 1, 2);
            return licznik / mianownik;
        } else {
            return Double.POSITIVE_INFINITY;
        }
    }

    double Skewness() {
        if (alpha > 3) {
            double licznik = (2 * r + alpha - 1) * (2 + beta + alpha - 1);
            double licznik2 = r * beta * (r + alpha - 1) * (beta + alpha - 1);
            double mianownik = (alpha - 3) * root(licznik2 / (alpha - 2), 2);
            return licznik / mianownik;
        } else {
            return Double.POSITIVE_INFINITY;
        }
    }
}
