
#include <bits.h> 
#include <iostream>
#include<cmath>
#include<map>
#include<vector>
#include<ctime>
#include<bits.h>
#include<unordered_set>

#include<math.h>

using namespace std;

long long int Fastexponential(long long int i, long long int j, long long int k) {

    long long int Base = i;
    long long int Power = j;
    long long int Mod = k;
    long long int Result = 1;

    long long int temp = 0;
    while (Power != 0)
    {
        if (Power % 2 == 0)
        {
            temp = Base * Base;
            Base = temp % Mod;
            Power = Power / 2;
            Result = Result;

        }
        if (Power % 2 == 1)
        {
            Result = (Result * Base) % Mod;
            Base = Base % Mod;
            Power = Power - 1;


        }

    }
    return Result;
}

bool miillerTest(long long int d, long long int n)
{
    long long int a = 2 + rand() % (n - 4);

    long long int x = Fastexponential(a, d, n);

    if (x == 1 || x == n - 1)
        return true;

    while (d != n - 1)
    {
        x = (x * x) % n;
        d *= 2;

        if (x == 1)      return false;
        if (x == n - 1)    return true;
    }

    return false;
}
 
bool isPrime(long long int n, long long int k)
{
    
    if (n <= 1 || n == 4)  return false;
    if (n <= 3) return true;

   
    long long int d = n - 1;
    while (d % 2 == 0)
        d /= 2;

    for (long long int i = 0; i < k; i++)
        if (!miillerTest(d, n))
            return false;

    return true;
}
long long int gcd(long long int a, long long int b)
{
	if (a == 0)
		return b;
	return gcd(b % a, a);
}
long long int gcdExtended(long long int a, long long int b, long long int* x, long long int* y)
{

	if (a == 0)
	{
		*x = 0;
		*y = 1;
		return b;
	}

	long long int x1, y1;

	long long int gcd = gcdExtended(b % a, a, &x1, &y1);


	*x = y1 - (b / a) * x1;
	*y = x1;

	if (*x < 0)
	{
		return *x + b;
	}
	return *x;

}

long long int BBGiant(long long int b, long long int t, long long int m) {
	long long int bound = floor(sqrt(m));
	long long int target = t;
	long long int base = b;
	long long int mod = m;
	map<  long long int, long long int> babaystep;
	map<  long long int, long long int>::iterator it;
	for (long long int i = 0; i < bound; i++) {
		babaystep.insert(pair< long long int, long long int>(i, Fastexponential(b, i, m)));

	}
	long long int x, y;
	long long int result = 0;
	long long int c = Fastexponential(gcdExtended(b, m, &x, &y), bound, m);
	long long int temp = t % m;
	for (long long int i = 1; i < bound; i++)
	{
		bool flag = false;
		temp = (temp * c) % m;

		for (it = babaystep.begin(); it != babaystep.end(); it++)
		{
			if (temp == it->second) {
				flag = true;

				result = (i * bound) + it->first;

				break;
			}

		}
		if (flag == true) {
			break;
		}
	}
	return result;
}


long long int PollardRho(long long int n) {
	long long int x, y, d;
	x = 2;
	y = 2;
	d = 1;
	while (d == 1)
	{
		x = (x * x + 1) % n;
		y = (y * y + 1) % n;
		y = (y * y + 1) % n;
		d = gcd(abs(x - y), n);
	}
	if (d == n) {
		return 0;
	}
	else {
		cout << d << "is a factor" << endl;
		return d;
	}
}

long long int Pollardp1(long long int p) {
	int bits = floor(log(p) / log(2));
	bool flag = false;
	while (flag == false) {
		srand(rand() % 100);
		long long int bound = sqrt(p);
		long long int b = rand()%bound- 1;
		if (gcd(b,p)==p)
		{
			Pollardp1(p);
		}
		if (gcd(b, p) == 1) {
			int p1 = 2;
			long long int temp = floor(log(p) / log(2));
			long long int L = 0;
			L = Fastexponential(b, pow(2, temp), p);
			long long int comp = gcd(L - 1, p);
				if (comp!=1)
				{
					std::cout << comp << "is a factor" << "such that"<<comp-1<<"is"<<b<<"-smooth"<<endl;
					return comp;
					
				}
				else
					for (long long int i = 3; i < p; i+=2)
					{
						if (isPrime(i, 6)) {
							temp = floor(log(p) / log(i));
							 L = Fastexponential(L, pow(i, temp), p);
							 comp = gcd(L - 1, p);
							if (comp != 1)
							{
								std::cout << comp << "is a factor" << "such that" << comp - 1 << "is" << b << "-smooth" << endl;
								return comp;

							}
						}
						else;
					}
			
		}
		
	}
	return 0;
}
long long int NaorR(long long int n) {
	long long int lowbound = pow(2, n);
	srand(unsigned(time(0)));
	long long int x = 1 + rand() % lowbound - 1;
	


	long long int p = 1;
	long long int q = 1;
	long long int count = 0;
	for (long long int i = 0; i < lowbound; i++) {
		if (count == 0) {
			p = rand() % lowbound + lowbound;
			if (isPrime(p, 6)) {
				count = 1;
			

			}
		}
		if (count == 1) {
			q = rand() % lowbound + lowbound;
			if (isPrime(q, 6)) {
				count = 2;
				

			}
		}
		if (count == 2)
		{
			break;
		}

	}
	long long int N = p * q;
	
	map< long long int, long long int>pairs;
	for (long long int i = 0; i < n + 1; i++)
	{
		long long int temp1 = rand() % N;
		long long int temp2 = rand() % N;
		pairs.insert(pair< long long int, long long int>(temp1, temp2));
	}
	vector< long long int> arr;
	for (long long int i = 0; i < n + 1; i++)
	{
		long long int temp = x;
		temp = temp >> i;
		temp = temp & 1;
		arr.push_back(temp);
	}
	
	map< long long int, long long int>::iterator itr;
	long long int sum = 0;
	long long int j = 0;
	for (itr = pairs.begin(); itr != pairs.end(); itr++) {
		if (arr[j] == 0) {
			sum += itr->second;
		}
		else {
			sum += itr->first;
		}
		j++;
	}
	
	bool ui = 0;
	long long int g = 0;
	while (!ui) {
		g = rand() % N;
		if (gcd(g, N) == 1) {
			ui = 1;
		}
	}
	long long int gin = 0;
	long long int temp4, temp5;
	gin = gcdExtended(g, N, &temp4, &temp5);

	long long int beta = Fastexponential(gin, sum, N);

	long long int r = 0;
	long long int rbound = pow(2, 2 * n + 1);
	r = rand() % rbound + rbound;
	unsigned long long int result;
	result = (r * beta) % 2;
	if (result==-1)
	{
		result + 2;
	}

	
	return result;

}

unordered_set <long long int> Primefactors1(long long int p) {
	unordered_set<long long int> factors;
	p = p - 1;
	long long int factor = 2;
	while (factor * factor <= p)
	{
		if (isPrime(factor, 6) && p % factor == 0)
		{
			factors.insert(factor);
			p = p / factor;
			factors.insert(p);
		}
		factor++;
	}
	
	
	return factors;
}


long long int Primitiverootsearch(unordered_set<long long int>s, long long int p) {
	bool flag = false;
	while (flag == false) {
		srand(rand() % 100);
		int j = rand() % 100;
		for (auto it = s.begin(); it != s.end(); it++) {
			long long int temp = (p - 1) / *it;

			long long int x = Fastexponential(j, temp, p);

			if (x == 1)
			{
				break;
			}

		}
		if (j > p)
		{
			j = j % p;
		}
		return j;
	}
	return 0;
}
bool Primitiverootverifier(long long int p, long long int r)
{
	bool flag = false;
		unordered_set<long long int> factors = Primefactors1(p);
		for (auto it = factors.begin(); it != factors.end(); it++) {
			long long int temp = (p - 1) / *it;

			long long int x = Fastexponential(r, temp, p);

			if (x == 1)
			{
				break;
				return false;
			}

		}
		if (r > p)
		{
			r = r % p;
		}
		flag = 1;
		return flag;
}
bool Primverifier(long long int p) {
	bool flag = false;
	flag = isPrime(p, 2);
	return flag;
}
long long int RamdomNumber(long long int bits) {
	long long int result = 0;
	vector< int>ranbits;
	srand(time(NULL));
	for (long long int f = 1; f < bits - 1; f++) {

		int temp = NaorR(bits - f);
		ranbits.push_back(temp);
	}
	result += pow(2, bits - 1);
	result += 1;
	for (int i = 0; i < ranbits.size(); i++)
	{
		result += pow(2, i + 1) * ranbits[i];
	}

	return result;
}
long long int RadomPrime(int bits) {
	bool flag = false;
	while (flag == false)
	{
		srand(rand() % 100);
		long long int j = RamdomNumber(bits);
		if (isPrime(j, 2 * bits) == true) {
			flag = 1;
			return j;
		}
	}
	return 0;
}


long long int ElgamalCipher(long long int text, long long int Pk, long long int p) {
	auto incrypted = text * Pk;
	incrypted = incrypted % p;
	if (incrypted < 0)
	{
		incrypted + p;
	}
	return incrypted;
}

void ELGAMALgenerator(int j) {
	long long int n = RadomPrime(j);
	long long int b = Primitiverootsearch(Primefactors1(n), n);
	cout << "n=" << n << endl;
	cout << "b=" << b << endl;

}
void ELGAMALAlice(long long int key, long long int base, long long int prime) {

	long long int AK = Fastexponential(base, key, prime);
	cout << "Alice's key:" << AK << endl;
}
void ELGAMGALBob(long long int key, long long int base, long long int prime) {
	bool flag1, flag2;
	flag1 = Primverifier(prime);
	flag2 = Primitiverootverifier(prime, base);
	if (flag1&&flag2)
	{
		long long int BK = Fastexponential(base, key, prime);
		cout << "Bob's key:" << BK << endl;
	}
	else
	{
		cout << "error!" << endl;
	}
}
void ELAliceencryption(long long int key, long long int base, long long int prime, long long int message) {
	long long int En = Fastexponential(base, key, prime);
	cout << "Encrypted Message is" << ElgamalCipher(message, En, prime) << endl;
}
void ELBobdecryption(long long int key, long long int base, long long int prime, long long int message) {
	long long int En = Fastexponential(base, key, prime);
	long long int x, y;
	long long int De = gcdExtended(En, prime, &x, &y);
	cout << "Decrypted Message is" << ElgamalCipher(message, De, prime) << endl;
}
void ELGAMALEVE(long long int Akey, long long int Bkey, long long int base, long long int prime, long long int message) {
	cout << "Im EVE" << endl;
	long long int A = BBGiant(base, Akey, prime);
	cout << "Alice's key is" << A << endl;
	long long int B = BBGiant(base, Bkey, prime);
	cout << "Bob's key is" << B << endl;
	long long int en = Fastexponential(Akey, B, prime);
	cout << "Encryption key is" << en << endl;
	long long int x, y;
	long long int DE = gcdExtended(en, prime, &x, &y);
	cout << "Decryption key is" << DE << endl;
	cout << "the message is" << ElgamalCipher(message, DE, prime);

}
void RSABobkeygenerate(int i,int j) {
	long long int p = RadomPrime(i);
	long long int q = RadomPrime(j);
	long long int n = p * q;
	long long int n1 = (p - 1) * (q - 1);
	bool flag = false;
	long long int e = 0;
	while (flag == false)
	{
		int i = 1;
		srand(rand() % 100);
		e = rand() % n;
		if (gcd(e, n1) == 1)
		{
			flag = 1;
		}
	}

	long long int x, y;
	long long int d = gcdExtended(e, n1, &x, &y);
	cout << "n=" << n << "e=" << e << "d=" << d << endl;


}
void RSAAliceencryption(long long int x, long long int e, long long int n) {
	auto en = Fastexponential(x, e, n);
	cout << "Encrypted Message is" << en << endl;
}
void RSABobdecryption(long long int x, long long int e, long long int n) {
	auto de = Fastexponential(x, e, n);
	cout << "Decrypted Message is" << de << endl;
}

void RSAEVE(long long int n,long long int e,long long int m) {
	cout << "I am EVE" << endl;
	long long int p = PollardRho(n);
	long long int q = n / p;
	cout << "p=" << p << endl << "q=" << q << endl;

	long long int n1 = (p - 1)*( q - 1);
	long long int x, y;
	long long int d = gcdExtended(e, n1, &x, &y);
	cout << "decryption key is" << d << endl;
	
		RSABobdecryption(m, d, n);

}

int main()
{
	
	cout << Primitiverootverifier(79487363, 11);
	
	return 0;

}