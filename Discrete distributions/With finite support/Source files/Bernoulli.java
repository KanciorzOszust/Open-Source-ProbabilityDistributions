import java.util.Random;

public class Bernoulli extends MathLibrary {
    double p;
    double q;
    Bernoulli(double p) {
        if (p < 0 || p > 1) {
            throw new IllegalArgumentException("0 <= p <= 1");
        }
        this.p = p;
        this.q = 1 - p;
    }

    double PMF(int k) {
        if (k == 0) return this.q;
        else if (k == 1) return this.p;
        throw new IllegalArgumentException("k = 0, k = 1");
    }

    double CDF(double k) {
        if (k < 0) return 0;
        else if (k > 0 && k < 1) return this.q;
        else return 1;
    }

    double Mean() {
        return this.p;
    }

    int Mediana() {
        if (p < 0.5) return 0;
        else if (p > 0.5) return 1;
        else {
            Random random = new Random();
            return random.nextInt(1);
        }
    }

    int Mode() {
        if (p < 0.5) return 0;
        else if (p > 0.5) return 1;
        else {
            Random random = new Random();
            return random.nextInt(1);
        }
    }

    double Variance() {
        return this.p * this.q;
    }

    double MAD() {
        return 0.5;
    }

    double Skewness() {
        return (this.q - this.p) / root(this.p * this.q, 2);
    }

    double EXkurtosis() {
        return (1 - (6 * this.p * this.q)) / (this.p * this.q);
    }

    double Entropy() {
        return (-1 * this.q * ln(this.q)) - (this.p * ln(this.p));
    }

    double MGF(double t) {
        return this.q + (this.p * euler(t));
    }

    double PGF(double z) {
        return q + (p * z);
    }

    double FisherInformation() {
        return 1.0 / (this.p * this.q);
    }
}
