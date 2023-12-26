public class Type1Gumbel extends MathLibrary {
    double mu;
    double beta;

    Type1Gumbel(double mu, double beta) {
        if (beta <= 0) throw new IllegalArgumentException("beta > 0");
        this.mu = mu;
        this.beta = beta;
    }

    private double z(double x) {
        return (x - mu) / beta;
    }

    double PDF(double x) {
        return (1 / beta) * euler(-1 * (z(x) + euler(-1 * z(x))));
    }

    double CDF(double x) {
        return euler(-1 * euler(-1 * z(x)));
    }

    double Mean() {
        return mu + (beta * EulerMascheroni);
    }

    double Median() {
        return mu - (beta * ln(2));
    }
    
    double Mode() {
        return mu;
    }

    double Variance() {
        return power(PI, 2) / 6 * power(beta, 2);
    }

    double Skewness() {
        return 12 * root(6, 2) * zeta(3) / power(PI, 3);
    }

    double EXkurtosis() {
        return 12.0 / 5;
    }

    double Entropy() {
        return ln(beta) + EulerMascheroni + 1;
    }

    double MGF(double t) {
        return gamma(1 - (beta * t)) * euler(mu * t);
    }
}