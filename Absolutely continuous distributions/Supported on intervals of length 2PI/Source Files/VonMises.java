public class VonMises extends MathLibrary{
    double mu;
    double k;

    VonMises(double mu, double k) {
        if (k <= 0) throw new IllegalArgumentException("k > 0");
        this.mu = mu;
        this.k = k;
    }

    double PDF(double x) {
        if (x < -1 * PI || x > PI) throw new IllegalArgumentException("-PI <= x <= PI");
        return euler(k * cosine(x - mu)) / (2 * PI * modifiedBessel(0, k));
    }

    double CDF(double x) {
        if (x < -1 * PI || x > PI) throw new IllegalArgumentException("-PI <= x <= PI");
        double value = 0;
        for (int j = 1; j < 7; j++) {
            value += modifiedBessel(j, k) * cosine(j * (x - mu));
        }
        return (1 / (2 * PI)) * (1 + ((2 / modifiedBessel(0, k)) * value));
    }

    double Mean() {
        return mu;
    }

    double Median() {
        return mu;
    }

    double Mode() {
        return mu;
    }

    double Variance() {
        return 1 - (modifiedBessel(1, k) / modifiedBessel(0, k));
    }

    double Entropy() {
        return (-1 * k * (modifiedBessel(1, k) / modifiedBessel(0, k))) + ln(2 * PI * modifiedBessel(0, k));
    }
}
