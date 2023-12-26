public class WrappedExponential extends MathLibrary{
    double lambda;

    WrappedExponential(double lambda) {
        if (lambda <= 0) throw new IllegalArgumentException("lambda > 0");
        this.lambda = lambda;
    }

    double PDF(double x) {
        if (x < 0 || x >= 2 * PI) throw new IllegalArgumentException("0 <= x < 2PI");
        double part1 = lambda * euler(-1 * lambda * x);
        double part2 = 1 - euler(-2 * PI * lambda);   
        return part1 / part2;
    }

    double CDF(double x) {
        if (x < 0 || x >= 2 * PI) throw new IllegalArgumentException("0 <= x < 2PI");
        double part1 = 1 - euler(-1 * lambda * x);
        double part2 = 1 - euler(-2 * PI * lambda);
        return part1 / part2;
    }

    double Mean() {
        return arcTangent(1 / lambda);
    }

    double Variance() {
        return 1 - (lambda / root(1 + power(lambda, 2), 2));
    }

    double Entropy() {
        double beta = euler(2 * PI * lambda);
        double part1 = 1 + ln((beta - 1) / lambda);
        double part2 = (beta / (beta - 1)) * ln(beta);
        return part1 - part2;
    }
}
