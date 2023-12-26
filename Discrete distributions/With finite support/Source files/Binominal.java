import java.util.Random;

public class Binominal extends MathLibrary {
    int n;
    double p;
    double q;
    Binominal(int n, double p) {
        if (p < 0 || p > 1) {
            throw new IllegalArgumentException("0 <= p <= 1");
        }
        this.n = n;
        this.p = p;
        this.q = 1 - n;
    }

    double PMF(int k) {
        return binomial(n, k) * power(p, k) * power(q, n - k);
    }

    double CDF(int k) {
        return regularizedIncompleteBeta(n - k, 1 + k, q);
    }

    double Mean() {
        return n * p;
    }

    int Median() {
        Random random = new Random();
        if (random.nextInt(1) == 0) return (int) (n * p);
        else return (int) (n * p + 1);
    }

    int Mode() {
        Random random = new Random();
        if (random.nextInt(1) == 0) {
            return (int) ((n + 1) * p);
        } else {
            return (int) ((n + 1) * p + 1) - 1;
        }
    }

    double Variance() {
        return n * p * q;
    }

    double Skewness() {
        return (q - p) / root(n * p * q, 2);
    }

    double EXkurtosis() {
        return (1 - (6 * p * q)) / (n * p * q);
    }

    double Entropy() {
        return 0.5 * log(17.07 * n * p * q, 2);
    }

    double MGF(double t) {
        return power(q + (p * euler(t)), n);
    }

    double PGF(double z) {
        return power(q + (p * z), n);
    }

    double FIsherInformation() {
        return n / (p * q);
    }
}