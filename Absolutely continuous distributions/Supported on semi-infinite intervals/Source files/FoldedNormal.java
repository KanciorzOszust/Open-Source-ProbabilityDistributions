public class FoldedNormal extends MathLibrary{
    double mu;
    double sigma;

    FoldedNormal(double mu, double sigma) {
        if (sigma <= 0) throw new IllegalArgumentException("sigma > 0");
        this.mu = mu;
        this.sigma = sigma;
    }

    private double NCDF(double x) {
        return 0.5 * (1 + errorFunction((x - mu) / root(2 * sigma, 2)));
    }

    double PDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        double part1 = (1 / root(2 * PI * sigma, 2)) * euler(-1 * power(x - mu, 2) / (2 * sigma));
        double part2 = (1 / root(2 * PI * sigma, 2)) * euler(-1 * power(x + mu, 2) / (2 * sigma));
        return part1 + part2;
    }

    double CDF(double x) {
        if (x < 0) throw new IllegalArgumentException("x >= 0");
        return 0.5 * (errorFunction((x + mu) / root(2 * sigma, 2)) + errorFunction((x - mu) / root(2 * sigma, 2)));
    }

    double Mean() {
        double part1 = root(sigma * (2 / PI), 2) * euler(power(-1 * mu, 2) / (2 * sigma));
        double part2 = mu * (1 - (2 * NCDF(-1 * mu / root(sigma, 2))));
        return part1 + part2;
    }
    
    double Variance() {
        return power(mu, 2) + sigma - power(Mean(), 2);
    }
}
