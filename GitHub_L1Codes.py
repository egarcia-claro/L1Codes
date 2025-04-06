 
    
    
#Algorithm 1 

def eta(r, m):
    if r==0:
        return 1
    
    n=len(m); indices = range(1, n+1); eta = 0
       
    for delta in range(r + 1):
        for J in Subsets(indices):
            sum_m_J = sum(m[i-1] for i in J)
            if delta - sum_m_J >= 0:
                binom_term = binomial(n + delta - sum_m_J - 1, delta - sum_m_J)
                eta += (-1) ** len(J) * binom_term
    
    return eta
        
 
 
   
#Algorithm 2 
    
def p_delta_o(delta, m):
    if delta==0:
        return 1
    
    n=len(m); l = [ceil((mi - 1) / 2) for mi in m]; indices=range(1, n+1);
    E = {i for i in indices if m[i-1] % 2 == 0};
    partial=sum(li for li in l); 
    p_delta = 0
    
    
    
    if len(E) != 0:
        raise ValueError("The m_i's must be odd")
    
    elif delta == partial:
        return 2**n
    
    P_n= Subsets(indices)
    for J in P_n:
        if len(J) > 0:
            for A in Subsets(J):
                sum_l_A = sum(l[i-1] + 1 for i in A)
                if delta - sum_l_A >= 0:
                    binom_term = binomial(len(J) + delta - sum_l_A - 1, delta - sum_l_A)
                    p_delta += (-1)**(n - len(J) + len(A)) * 2**len(J) * binom_term

    return p_delta




#Algorithm 3 

def p_delta_e(delta, m):
    if delta==0:
        return 1
    
    n=len(m); l = [ceil((mi - 1) / 2) for mi in m]; indices=range(1, n+1);
    E = {i for i in indices if m[i-1] % 2 == 0};
    partial=sum(li for li in l); 
    p_delta = 0
    
    
    
    if len(E) != n:
        raise ValueError("The m_i's must be even")

    if delta == partial:
        return 1
    
    P_n = Subsets(indices)
    for J in P_n:
        if len(J) > 0:  # Ensure J is non-empty
            J_c = set(indices) - set(J)
            u_J_delta = (-1)**(len(J)) if delta == sum(l[i-1] for i in J_c) else 0
            
            term_sum = u_J_delta
            
            for A in Subsets(J):
                if len(A) > 0:
                    for B in Subsets(A):
                        sum_l_B_J_c = sum(l[i-1] for i in B) + sum(l[i-1] for i in J_c)
                        if delta - sum_l_B_J_c >= 0:
                            binom_term = binomial(len(A) + delta - sum_l_B_J_c - 1, delta - sum_l_B_J_c)
                            term_sum += (-1)**(len(J) - len(A) + len(B))  *  2**len(A)  *  binom_term
            
            p_delta += term_sum
    
    return p_delta




#Algorithm 4 

def p_delta(delta, m):
    if delta==0:
        return 1
    
    n=len(m); indices=range(1, n+1)
    m_e=[m[i-1] for i in indices if m[i-1] % 2 == 0]
    partial_e=sum(ceil((m - 1) / 2) for m in m_e)
    m_o=[m[i-1] for i in indices if m[i-1] % 2 != 0]
    p_delta = 0
    
    return sum(p_delta_e(j, m_e) *  p_delta_o(delta-j, m_o) for j in range(partial_e + 1))
    



#Algorithm 5

def gamma(r, m):
    if r==0:
        return 1
    
    return 1 + sum(p_delta(delta, m) for delta in range(1, r+1))





#Algorithm 6

def sphere_size(j, x, m):
    l_min = max(x, m - (x + 1))
    l_max = min(x, m - (x + 1))
    
    if j == 0:
        return 1
    elif 0 < j <= l_max:
        return 2
    elif l_max < j <= l_min:
        return 1
    else:
        return 0





#Algorithm 7

def compute_p_t(t,a, m_list):
    product = 1
    for i in range(len(a)):
        m_i = m_list[i]
        l_i = max(a[i], m_i - (a[i] + 1))
        sum_expr = sum(sphere_size(j, a[i], m_i) * t**j for j in range(l_i + 1))
        product *= sum_expr
    return product





#Algorithm 8

def generate_representative_orbits(m,n):
    from itertools import product

    bounds = list(range(0, (m - 1) // 2 + 1)); reps = [ ] 
    cartesian_product=list(product(bounds, repeat=n))

    for candidate in cartesian_product:
        if list(candidate) == sorted(candidate):
            reps.append(candidate)
    return reps





#Algorithm 9

def compute_stabilizer_size(a, m):
       
    from collections import Counter
    from math import factorial

    n = len(a);    I_a = set()
    for i in range(n):
        reflected_value = m - (a[i] + 1)
        if reflected_value in a:
            I_a.add(i)

    multiplicities = Counter(a)  # Dictionary with value: count
    product_factorial = 1
    for count in multiplicities.values():
        product_factorial *= factorial(count)

    stabilizer_size = (2 ** len(I_a)) * product_factorial
    return stabilizer_size





#Algorithm 10

def compute_p_bar(t, m, n):
    from math import factorial
    T = generate_representative_orbits(m,n); total_sum = 0;
    for a in T:
        stab_size = compute_stabilizer_size(a, m)
        orbit_size = QQ(factorial(len(a)) * (2**len(a))) / stab_size
        p_t_a = compute_p_t(t, a, [m]*n)
        total_sum += orbit_size * p_t_a

    return total_sum / m**n




    
    






  