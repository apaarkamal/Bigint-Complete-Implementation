// wherever you include this header file
// convert you int main to int32_t main()
// because i defined int as long long int

// example
// #include "bigint.h"
// int32_t main(){
//  bigint a,b;
//	cin>>a>>b;
//	cout<<a*b;
// }

#define int long long int

const int base = 1000000000;
const int base_digits = 9;

struct bigint {
	// this vector carries at most 9 digit numbers at one index
	vector<int> z;

	// sign of the number
	//  1  positive
	// -1  negative
	int sign;

	bigint() : sign(1) {}

	bigint(int v) {
		*this = v;
	}

	bigint(const string &s) {
		read(s);
	}

	void operator=(const bigint &v) {
		sign = v.sign;
		z = v.z;
	}

	void operator=(int v) {
		sign = 1;
		if (v < 0) {
			sign = -1, v = -v;
		}
		z.clear();
		for (; v > 0; v = v / base) {
			z.push_back(v % base);
		}
	}

	bigint operator+(const bigint &v) const {
		if (sign == v.sign) {
			bigint res = v;
			for (int i = 0, carry = 0; i < (int) max(z.size(), v.z.size()) || carry; ++i) {
				if (i == (int) res.z.size()) {
					res.z.push_back(0);
				}
				res.z[i] += carry + (i < (int) z.size() ? z[i] : 0);
				carry = res.z[i] >= base;
				if (carry) {
					res.z[i] -= base;
				}
			}
			return res;
		}
		else {
			return *this - (-v);
		}
	}

	bigint operator-(const bigint &v) const {
		if (sign == v.sign) {
			if (abs() >= v.abs()) {
				bigint res = *this;
				for (int i = 0, carry = 0; i < (int) v.z.size() || carry; ++i) {
					res.z[i] -= carry + (i < (int) v.z.size() ? v.z[i] : 0);
					carry = res.z[i] < 0;
					if (carry) {
						res.z[i] += base;
					}
				}
				res.trim();
				return res;
			}
			else {
				return -(v - *this);
			}
		}
		else {
			return *this + (-v);
		}
	}

	void operator*=(int v) {
		if (v < 0) {
			sign = -sign, v = -v;
		}
		for (int i = 0, carry = 0; i < (int) z.size() || carry; ++i) {
			if (i == (int) z.size()) {
				z.push_back(0);
			}
			int cur = z[i] * v + carry;
			carry = cur / base;
			z[i] = cur % base;
		}
		trim();
	}

	bigint operator*(int v) const {
		bigint res = *this;
		res *= v;
		return res;
	}

	friend pair<bigint, bigint> divmod(const bigint &a1, const bigint &b1) {
		int norm = base / (b1.z.back() + 1);
		bigint a = a1.abs() * norm;
		bigint b = b1.abs() * norm;
		bigint q, r;
		q.z.resize(a.z.size());
		for (int i = a.z.size() - 1; i >= 0; i--) {
			r *= base;
			r += a.z[i];
			int s1 = b.z.size() < r.z.size() ? r.z[b.z.size()] : 0;
			int s2 = b.z.size() - 1 < r.z.size() ? r.z[b.z.size() - 1] : 0;
			int d = (s1 * base + s2) / b.z.back();
			r -= b * d;
			while (r < 0) {
				r += b, --d;
			}
			q.z[i] = d;
		}
		q.sign = a1.sign * b1.sign;
		r.sign = a1.sign;
		q.trim();
		r.trim();
		return make_pair(q, r / norm);
	}

	bigint operator/(const bigint &v) const {
		return divmod(*this, v).first;
	}

	bigint operator%(const bigint &v) const {
		return divmod(*this, v).second;
	}

	void operator/=(int v) {
		if (v < 0) {
			sign = -sign, v = -v;
		}
		for (int i = z.size() - 1, rem = 0; i >= 0; --i) {
			int cur = z[i] + rem * base;
			z[i] = cur / v;
			rem = cur % v;
		}
		trim();
	}

	bigint operator/(int v) const {
		bigint res = *this;
		res /= v;
		return res;
	}

	int operator%(int v) const {
		if (v < 0) {
			v = -v;
		}
		int m = 0;
		for (int i = z.size() - 1; i >= 0; --i) {
			m = (m * base + z[i]) % v;
		}
		return m * sign;
	}

	void operator+=(const bigint &v) {
		*this = *this + v;
	}
	void operator-=(const bigint &v) {
		*this = *this - v;
	}
	void operator*=(const bigint &v) {
		*this = *this * v;
	}
	void operator/=(const bigint &v) {
		*this = *this / v;
	}
	void operator%=(const bigint &v) {
		*this = *this % v;
	}
	void operator%=(int v) {
		*this = *this % v;
	}

	bool operator<(const bigint &v) const {
		if (sign != v.sign) {
			return sign < v.sign;
		}
		if (z.size() != v.z.size()) {
			return z.size() * sign < v.z.size() * v.sign;
		}
		for (int i = z.size() - 1; i >= 0; i--) {
			if (z[i] != v.z[i]) {
				return z[i] * sign < v.z[i] * sign;
			}
		}
		return false;
	}

	bool operator>(const bigint &v) const {
		return v < *this;
	}
	bool operator<=(const bigint &v) const {
		return !(v < *this);
	}
	bool operator>=(const bigint &v) const {
		return !(*this < v);
	}
	bool operator==(const bigint &v) const {
		return !(*this < v) && !(v < *this);
	}
	bool operator!=(const bigint &v) const {
		return *this < v || v < *this;
	}

	void trim() {
		while (!z.empty() && z.back() == 0) {
			z.pop_back();
		}
		if (z.empty()) {
			sign = 1;
		}
	}

	bool isZero() const {
		return z.empty() || ((int) z.size() == 1 && !z[0]);
	}

	bigint operator-() const {
		bigint res = *this;
		res.sign = -sign;
		return res;
	}

	bigint abs() const {
		bigint res = *this;
		res.sign *= res.sign;
		return res;
	}

	int longValue() const {
		int res = 0;
		for (int i = z.size() - 1; i >= 0; i--) {
			res = res * base + z[i];
		}
		return res * sign;
	}

	friend bigint gcd(const bigint &a, const bigint &b) {
		return b.isZero() ? a : gcd(b, a % b);
	}
	friend bigint lcm(const bigint &a, const bigint &b) {
		return a / gcd(a, b) * b;
	}

	void read(const string &s) {
		sign = 1;
		z.clear();
		int pos = 0;
		while (pos < (int) s.size() && (s[pos] == '-' || s[pos] == '+')) {
			if (s[pos] == '-') {
				sign = -sign;
			}
			++pos;
		}
		for (int i = s.size() - 1; i >= pos; i -= base_digits) {
			int x = 0;
			for (int j = max(pos, i - base_digits + 1); j <= i; j++) {
				x = x * 10 + s[j] - '0';
			}
			z.push_back(x);
		}
		trim();
	}

	friend istream &operator>>(istream &stream, bigint &v) {
		string s;
		stream >> s;
		v.read(s);
		return stream;
	}

	friend ostream &operator<<(ostream &stream, const bigint &v) {
		if (v.sign == -1) {
			stream << '-';
		}
		stream << (v.z.empty() ? 0 : v.z.back());
		for (int i = v.z.size() - 2; i >= 0; --i) {
			stream << setw(base_digits) << setfill('0') << v.z[i];
		}
		return stream;
	}

	static vector<int> convert_base(const vector<int> &a, int old_digits, int new_digits) {
		vector<int> p(max(old_digits, new_digits) + 1);
		p[0] = 1;
		for (int i = 1; i < (int) p.size(); i++) {
			p[i] = p[i - 1] * 10;
		}
		vector<int> res;
		int cur = 0;
		int cur_digits = 0;
		for (int i = 0; i < (int) a.size(); i++) {
			cur += a[i] * p[cur_digits];
			cur_digits += old_digits;
			while (cur_digits >= new_digits) {
				res.push_back(cur % p[new_digits]);
				cur /= p[new_digits];
				cur_digits -= new_digits;
			}
		}
		res.push_back(cur);
		while (!res.empty() && res.back() == 0) {
			res.pop_back();
		}
		return res;
	}

	typedef vector<int> vll;

	static vll karatsubaMultiply(const vll &a, const vll &b) {
		int n = a.size();
		vll res(n + n);
		if (n <= 32) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					res[i + j] += a[i] * b[j];
				}
			}
			return res;
		}
		int k = n >> 1;
		vll a1(a.begin(), a.begin() + k);
		vll a2(a.begin() + k, a.end());
		vll b1(b.begin(), b.begin() + k);
		vll b2(b.begin() + k, b.end());
		vll a1b1 = karatsubaMultiply(a1, b1);
		vll a2b2 = karatsubaMultiply(a2, b2);
		for (int i = 0; i < k; i++) {
			a2[i] += a1[i];
		}
		for (int i = 0; i < k; i++) {
			b2[i] += b1[i];
		}
		vll r = karatsubaMultiply(a2, b2);
		for (int i = 0; i < (int) a1b1.size(); i++) {
			r[i] -= a1b1[i];
		}
		for (int i = 0; i < (int) a2b2.size(); i++) {
			r[i] -= a2b2[i];
		}
		for (int i = 0; i < (int) r.size(); i++) {
			res[i + k] += r[i];
		}
		for (int i = 0; i < (int) a1b1.size(); i++) {
			res[i] += a1b1[i];
		}
		for (int i = 0; i < (int) a2b2.size(); i++) {
			res[i + n] += a2b2[i];
		}
		return res;
	}

	bigint operator*(const bigint &v) const {
		vector<int> a6 = convert_base(this->z, base_digits, 6);
		vector<int> b6 = convert_base(v.z, base_digits, 6);
		vll a(a6.begin(), a6.end());
		vll b(b6.begin(), b6.end());
		while (a.size() < b.size()) {
			a.push_back(0);
		}
		while (b.size() < a.size()) {
			b.push_back(0);
		}
		while (a.size() & (a.size() - 1)) {
			a.push_back(0);
			b.push_back(0);
		}
		vll c = karatsubaMultiply(a, b);
		bigint res;
		res.sign = sign * v.sign;
		for (int i = 0, carry = 0; i < (int) c.size(); i++) {
			int cur = c[i] + carry;
			res.z.push_back(cur % 1000000);
			carry = cur / 1000000;
		}
		res.z = convert_base(res.z, 6, base_digits);
		res.trim();
		return res;
	}
};

template <typename T>
T power(T a, long long b) {
	T r = 1;
	while (b) {
		if (b & 1) {
			r *= a;
		}
		a *= a;
		b >>= 1;
	}
	return r;
}

template <typename T>
T inverse(T a, T m) {
	a = a % m;
	if (a < 0) {
		a += m;
	}
	T b = m, u = 0, v = 1;
	while (a != 0) {
		T t = b / a;
		b -= a * t;
		swap(a, b);
		u -= v * t;
		swap(u, v);
	}
	if (u < 0) {
		u += m;
	}
	return u;
}
