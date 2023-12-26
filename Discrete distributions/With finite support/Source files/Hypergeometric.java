public class Hypergeometric extends MathLibrary{
    int N;
    int K;
    int n;
    
    Hypergeometric(int N, int K, int n) {
        if (N < 0 || K < 0 || n < 0) throw new IllegalArgumentException("N >= 0, K >= 0, n >= 0");
        if (K > N || n > N) throw new IllegalArgumentException("K < N, n < N");
        this.N = N;
        this.K = K;
        this.n = n;
    }

    double PMF(int k) {
        if (k < max(0, n + K - N) || k > min(n, K)) throw new IllegalArgumentException("max(0, n + K - N) <= k <= min(n , K)");
        return (binomial(K, k) * binomial(N - K, n - k)) / binomial(N, n);
    }

    double Mean() {
        return n * ((double) K / N);
    }

    int Mode() {
        return (int) (((n + 1) * (K + 1)) / (N + 2));
    }

    double Variance() {
        return n * (((double) K / N)) * ((double) (N - K) / N) * ((double) (N - n) / N - 1); 
    }

    double Skewness() {
        double part1 = (N - (2 * K)) * root(N - 1, 2) * (N - (2 * n));
        double part2 = root(n * K * (N - K) * (N - n), 2) * (N - 2);
        return part1 / part2;
    }

    double EXkurtosis() {
        double part1 = 1.0 / (n * K * (N - K) * (N - n) * (N - 2) * (N - 3));
        double part2 = (N - 1) * power(N, 2) * (N * (N + 1) - (6 * K * (N - K)) - (6 * n * (N - n)));
        double part3 = 6 * n * K * (N - K) * (N - n) * (5 * N - 6);
        return part1 * (part2 + part3);
    }

    double MGF(double t) {
        double part = binomial(N - K, n) * hypergeometric(-1 * n, -1 * K, N - K - n, euler(t));
        return part / binomial(N, n);
    }
}
