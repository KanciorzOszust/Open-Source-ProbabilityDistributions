public class InverseGaussian extends MathLibrary{
    double mu;
    double gamma;

    InverseGaussian(double mu, double gamma) {
        if (mu <= 0 || gamma <= 0) throw new IllegalArgumentException("mu, gamma > 0");
        this.mu = mu;
        this.gamma = gamma;
    }

    double PDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        double part1 = root(gamma / (2 * PI * power(x, 3)), 2);
        double part2 = euler(-1 * (gamma * power(x - mu, 2)) / (2 * x * power(mu, 2)));
        return part1 * part2;
    }

    private double PHI(double x) {
        return 0.5 * (1 + errorFunction(x / root(2, 2)));
    }

    double CDF(double x) {
        if (x <= 0) throw new IllegalArgumentException("x > 0");
        double part1 = PHI(power(gamma / x, 2) * (x / mu - 1));
        double part2 = PHI(-1 * root(gamma / x, 2) * (x / mu + 1));
        return part1 + (euler((2 * gamma) / mu) * part2); 
    }

    double Mean() {
        return mu;
    }

    double Mode() {
        return mu * (root(1 + ((9 * power(mu, 2)) / (4 * power(gamma, 2))), 2) - ((3 * mu) / (2 * gamma)));
    }

    double Variance() {
        return power(mu, 3) / gamma;
    }

    double Skewness() {
        return 3 * root(mu / gamma, 2);
    }

    double EXkurtosis() {
        return (15 * mu) / gamma;
    }

    double MGF(double t) {
        double part = 1 - root(1 - ((2 * t * power(mu, 2)) / gamma), 2);
        return euler((gamma / mu) * part);
    }
}
