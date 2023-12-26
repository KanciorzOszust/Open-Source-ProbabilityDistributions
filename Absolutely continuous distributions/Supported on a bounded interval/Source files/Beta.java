import java.util.Random;

public class Beta extends MathLibrary{
    double alpha;
    double beta;

    Beta(double alpha, double beta) {
        if (alpha <= 0 || beta <= 0) throw new IllegalArgumentException("alpha, beta > 0");
        this.alpha = alpha;
        this.beta = beta;
    }

    double PDF(double x) { 
        return (doublePower(x, alpha - 1) * doublePower(1 - x, beta - 1)) / beta(alpha, beta);
    }

    double CDF(double x) {
        return regularizedIncompleteBeta(alpha, beta, x);
    }

    double Mean() {
        return alpha / (alpha + beta);
    }

    double Median() {
        return 1.0 / regularizedIncompleteBeta(alpha, beta, 0.5);
    }

    double Mode() {
        Random random = new Random();
        if (alpha > 1 && beta > 1) {
            return (alpha - 1) / (alpha + beta - 2);
        } else if (alpha == 1 && beta == 1) {
            return random.nextDouble();
        } else if (alpha < 1 && beta < 1) {
            return random.nextInt(1);
        } else if (alpha <= 1 && beta > 1) {
            return 0;
        } else {
            return 1;
        }
    }

    double Variance() {
        return (alpha * beta) / (power(alpha + beta, 2) * (alpha + beta + 1));
    }

    double Skewness() {
        return (2 * (beta - alpha) * root(alpha + beta + 1, 2)) / ((alpha + beta +2) * root(alpha * beta, 2));
    }

    double EXkurtosis() {
        double licznik = 6 * (power(alpha - beta, 2) * (alpha + beta + 1) - (alpha * beta * (alpha + beta + 2)));
        double mianownik = alpha * beta * (alpha + beta + 2) * (alpha + beta + 3);
        return licznik / mianownik;
    }

    double Entropy() {
        return ln(beta(alpha, beta)) - ((alpha - beta) * digamma(alpha)) - ((beta - 1) * digamma(beta)) + ((alpha + beta - 2) * digamma(alpha + beta));
    }

    private double MGFpart1(int k) {
        double wartosc = 1;
        for (int r = 0; r < k - 1; r++) { 
            wartosc *= (alpha + r) / (alpha + beta + r);
        }
        return wartosc;
    }

    double MGF(double t) {
        double wartosc = 0;
        for (int k = 1; k < 7; k++) {
            wartosc += MGFpart1(k) * power(t, k) * factorial(k);
        }
        return 1 + wartosc;
    }
}
