public class Skellam extends MathLibrary{
    double mu1;
    double mu2;

    Skellam(double mu1, double mu2) {
        if (mu1 < 0 || mu2 < 0) throw new IllegalArgumentException("mu1 >= 0, mu2 >= 0");
        this.mu1 = mu1;
        this.mu2 = mu2;
    }

    double Mean() {
        return mu1 - mu2;
    }

    double Variance() {
        return mu1 + mu2;
    }

    double Skewness() {
        return (mu1 - mu2) / ((mu1 + mu2) * root(mu1 + mu2, 2));
    }

    double EXkurtosis() {
        return 1 / (mu1 + mu2);
    }

    double MGF(double t) {
        return euler(-1 * (mu1 + mu2) + (mu1 * euler(t)) + (mu2 * euler(-1 * t)));
    }
}
