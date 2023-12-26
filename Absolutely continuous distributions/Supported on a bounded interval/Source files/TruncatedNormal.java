public class TruncatedNormal extends MathLibrary{
    double mu;
    double sigma;
    double a;
    double b;
    double eta;
    double alpha;
    double beta;
    double Z;

    TruncatedNormal(double mu, double sigma, double a, double b) {
        if (sigma < 0) throw new IllegalArgumentException("sigma >= 0");
        if (sigma <= a && sigma >= b) throw new IllegalArgumentException("a < sigma < b");
        if (b <= a) throw new IllegalArgumentException("b > a");
        this.mu = mu;
        this.sigma = root(sigma, 2);
        this.a = a;
        this.b = b;
        this.alpha = (a - mu) / root(sigma, 2);
        this.beta = (b - mu) / root(sigma, 2);
        this.Z = PHI(beta) - PHI(alpha);
    }

    private double eta(double x) {
        return (x - mu) / sigma;
    }

    private double phi(double x) {
        return (1 / root(2 * PI, 2)) * euler(-0.5 * power(x, 2));
    }

    private double PHI(double x) {
        return 0.5 * (1 + errorFunction(x / root(2, 2)));
    }

    double PDF(double x) {
        if (x < a || x > b) throw new IllegalArgumentException("a <= x <= b");
        return phi(eta(x)) / (sigma * Z);
    }

    double CDF(double x) {
        if (x < a || x > b) throw new IllegalArgumentException("a <= x <= b");
        return (PHI(eta(x)) - PHI(alpha)) / Z;
    }

    double Mean() {
        return mu + (((phi(alpha) - phi(beta)) / Z) * sigma);
    }

    double Mode() {
        if (mu < a) return a;
        else if (a <= mu && mu <= b) return mu;
        else return b;
    }

    double Variance() {
        double part1 = 1 - ((beta * phi(beta) - (alpha * phi(alpha))) / Z);
        double part2= power((phi(alpha) - phi(beta)) / Z, 2);
        return power(sigma, 2) * (part1 - part2);
    }

    double Entropy() {
        double part1 = ln(root(2 * PI * euler(1), 2) * sigma * Z);
        double part2 = ((alpha * phi(alpha)) - (beta * phi(beta))) / (2 * Z);
        return part1 + part2;
    }

    double MGF(double t) {
        double part1 = euler((mu * t + (power(sigma, 2) * power(t, 2))) / 2);
        double part2 = (PHI(beta - (sigma * t)) - PHI(alpha - (sigma * t))) / (PHI(beta) - PHI(alpha));
        return part1 * part2;
    }
}
