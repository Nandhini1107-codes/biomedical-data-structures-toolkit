"""
U25ADG02 - Data Structures and Algorithms | Group 7 Project
Module 1 (Unit 3) T1 : Expression tree for clinical formulas (BMI, BSA, MAP)
Module 2 (Unit 4) T8 : BFS for contact tracing
Module 3 (Unit 5) T15 : Recursive vs iterative bacterial growth model
All data is synthetic. Core structures (stack, queue, tree, graph) are written by us.
Run: python G07_project.py
"""
import time
import sys

sys.setrecursionlimit(5000)


# =====================================================================
# Basic helper structures (own implementations, no ready-made containers)
# =====================================================================
class Stack:
    """Linked-list stack. push/pop/peek: O(1) time, O(n) space overall."""
    class _N:
        def __init__(self, v, nxt):
              self.v, self.next = v, nxt

    def __init__(self):
        self.top = None

    def push(self, v):
        self.top = Stack._N(v, self.top)

    def pop(self):
        if self.top is None:
            raise IndexError("stack empty")
        v = self.top.v
        self.top = self.top.next
        return v

    def peek(self):
        return None if self.top is None else self.top.v

    def is_empty(self):
        return self.top is None


class Queue:
    """Linked-list queue. enqueue/dequeue: O(1) time."""
    class _N:
        def __init__(self, v):
              self.v, self.next = v, None

    def __init__(self):
        self.head = self.tail = None

    def enqueue(self, v):
        n = Queue._N(v)
        if self.tail is None:
            self.head = self.tail = n


        else:
            self.tail.next = n
            self.tail = n

    def dequeue(self):
        if self.head is None:
            raise IndexError("queue empty")
        v = self.head.v
        self.head = self.head.next
        if self.head is None:
            self.tail = None
        return v

    def is_empty(self):
        return self.head is None


# =====================================================================
# MODULE 1 : Expression tree (Unit 3, T1)
# =====================================================================
PREC = {'+': 1, '-': 1, '*': 2, '/': 2, '^': 3}
RIGHT_ASSOC = {'^'}

# Synthetic clinical formula library (infix strings)
FORMULAS = {


    '1': ("BMI = weight(kg) / height(m)^2", "w / h ^ 2", ['w', 'h']),
    '2': ("BSA (Mosteller) = sqrt(h(cm)*w(kg)/3600)", "( h * w / 3600 ) ^ 0.5", ['h', 'w']),
    '3': ("MAP = (SBP + 2*DBP) / 3", "( s + 2 * d ) / 3", ['s', 'd']),
}


class TNode:
    def __init__(self, val, left=None, right=None):
        self.val, self.left, self.right = val, left, right


def tokenize(expr):
    """Split an infix string into numbers, variables, operators, brackets. O(n)."""
    tokens, i, n = [], 0, len(expr)
    while i < n:
        c = expr[i]
        if c == ' ':
            i += 1
        elif c.isdigit() or c == '.':
            j = i
            while j < n and (expr[j].isdigit() or expr[j] == '.'):
                  j += 1
            tokens.append(expr[i:j])
            i = j
        elif c.isalpha():
            j = i
            while j < n and (expr[j].isalnum() or expr[j] == '_'):
                  j += 1
            tokens.append(expr[i:j])
            i = j
        elif c in PREC or c in '()':
            tokens.append(c)
            i += 1
        else:
            raise ValueError("Invalid character '%s'" % c)
    return tokens


def infix_to_postfix(tokens):
    """Shunting-yard. Time O(n), space O(n)."""
    out, st = [], Stack()
    for t in tokens:
        if t in PREC:
            while (not st.is_empty() and st.peek() in PREC and
                   (PREC[st.peek()] > PREC[t] or
                      (PREC[st.peek()] == PREC[t] and t not in RIGHT_ASSOC))):
                out.append(st.pop())
            st.push(t)
        elif t == '(':
            st.push(t)
        elif t == ')':
            while not st.is_empty() and st.peek() != '(':
                out.append(st.pop())


            if st.is_empty():
                raise ValueError("Mismatched brackets")
            st.pop()
        else:
            out.append(t)
    while not st.is_empty():
        t = st.pop()
        if t == '(':
            raise ValueError("Mismatched brackets")
        out.append(t)
    return out


def build_tree(postfix):
    """Build expression tree from postfix. Time O(n), space O(n)."""
    st = Stack()
    for t in postfix:
        if t in PREC:
            if st.is_empty():
                 raise ValueError("Malformed expression")
            r = st.pop()
            if st.is_empty():
                 raise ValueError("Malformed expression")
            l = st.pop()
            st.push(TNode(t, l, r))
        else:
            st.push(TNode(t))
    if st.is_empty():
        raise ValueError("Empty expression")
    root = st.pop()
    if not st.is_empty():
        raise ValueError("Malformed expression")
    return root


def inorder(n):
    if n is None:
        return ''
    if n.left is None and n.right is None:
        return n.val
    return '(' + inorder(n.left) + ' ' + n.val + ' ' + inorder(n.right) + ')'


def preorder(n):
    if n is None:
        return ''
    return (n.val + ' ' + preorder(n.left) + ' ' + preorder(n.right)).strip()


def postorder(n):
    if n is None:
        return ''
    return (postorder(n.left) + ' ' + postorder(n.right) + ' ' + n.val).strip()


def evaluate(n, env):
    """Post-order evaluation. Time O(n), space O(h) recursion."""
    if n.left is None and n.right is None:
        try:
             return float(n.val)
        except ValueError:
             if n.val not in env:
                 raise KeyError("No value for variable '%s'" % n.val)
             return env[n.val]
    a, b = evaluate(n.left, env), evaluate(n.right, env)
    if n.val == '+': return a + b
    if n.val == '-': return a - b
    if n.val == '*': return a * b
    if n.val == '/':
        if b == 0:
             raise ZeroDivisionError("Division by zero")
        return a / b
    if n.val == '^':
        if a < 0 and b != int(b):
             raise ValueError("Negative base with fractional power")
        return a ** b


def height(n):
    return 0 if n is None else 1 + max(height(n.left), height(n.right))


def print_tree(n, level=0):
    if n is None:
        return
    print_tree(n.right, level + 1)
    print('       ' * level + n.val)
    print_tree(n.left, level + 1)


def read_float(prompt, lo=None, hi=None):
    while True:
        try:
             v = float(input(prompt))
             if (lo is not None and v < lo) or (hi is not None and v > hi):
                 print(" Value must be %s%s%s." % (
                      ">= %s" % lo if lo is not None else "",
                      " and " if lo is not None and hi is not None else "",
                      "<= %s" % hi if hi is not None else ""))
                 continue
             return v
        except ValueError:
             print(" Invalid number, try again.")


def module1():
    while True:
        print("\n--- MODULE 1: Expression Tree for Clinical Formulas ---")
        print("1. BMI   2. BSA   3. MAP   4. Custom infix expression   0. Back")
        ch = input("Choice: ").strip()
        if ch == '0':
            return
        if ch in FORMULAS:
            title, infix, vars_ = FORMULAS[ch]
            print("Formula:", title)
        elif ch == '4':
            infix = input("Enter infix (use spaces optional, ops + - * / ^): ")
            title = "custom"
        else:
            print("Invalid choice.")
            continue


        try:
            toks = tokenize(infix)
            pf = infix_to_postfix(toks)
            root = build_tree(pf)
        except ValueError as e:
            print("Error:", e)
            continue
        print("Postfix :", ' '.join(pf))
        print("Inorder :", inorder(root))
        print("Preorder :", preorder(root))
        print("Postorder:", postorder(root))
        print("Tree height:", height(root))
        print("Tree (rotated, root at left):")
        print_tree(root)
        names = []
        for t in toks:
            if t not in PREC and t not in '()' and not t.replace('.', '').isdigit() and t not in names:
                names.append(t)
        env = {}
        for v in names:
            env[v] = read_float("Enter value for %s: " % v, 0 if ch != '4' else None)
        try:
            print("Result = %.4f" % evaluate(root, env))
        except (ZeroDivisionError, KeyError, ValueError) as e:
            print("Error:", e)


# =====================================================================
# MODULE 2 : BFS contact tracing (Unit 4, T8)
# =====================================================================
class ContactGraph:
    """Undirected graph, adjacency list. Space O(V+E)."""

    def __init__(self):
        self.names = []      # index -> id
        self.adj = []        # index -> list of neighbour indices


    def index_of(self, name):
        for i in range(len(self.names)):
            if self.names[i] == name:
                return i
        return -1

    def add_person(self, name):
        if self.index_of(name) != -1:
            return False
        self.names.append(name)
        self.adj.append([])
        return True

    def add_contact(self, a, b):
        """O(V) because of lookup + duplicate check."""
        i, j = self.index_of(a), self.index_of(b)
        if i == -1 or j == -1 or i == j:
            return False
        if j in self.adj[i]:
            return False
        self.adj[i].append(j)
        self.adj[j].append(i)
        return True

    def bfs(self, src):
        """Returns (level, parent). Time O(V+E), space O(V)."""
        V = len(self.names)
        level = [-1] * V
        parent = [-1] * V
        q = Queue()
        level[src] = 0
        q.enqueue(src)
        while not q.is_empty():
            u = q.dequeue()
            for v in self.adj[u]:
                if level[v] == -1:
                    level[v] = level[u] + 1
                    parent[v] = u
                    q.enqueue(v)
        return level, parent

    def chain(self, a, b):
        """Shortest chain a -> b (fewest contacts) or None."""
        i, j = self.index_of(a), self.index_of(b)
        if i == -1 or j == -1:
            return None
        level, parent = self.bfs(i)
        if level[j] == -1:
            return None
        path, cur = [], j


        while cur != -1:
            path.append(self.names[cur])
            cur = parent[cur]
        return path[::-1]


def sample_graph():
    """24 synthetic persons (P01..P24); P22-P24 form a separate cluster."""
    g = ContactGraph()
    for k in range(1, 25):
        g.add_person("P%02d" % k)
    edges = [(1, 2), (1, 3), (1, 4), (2, 5), (2, 6), (3, 6), (3, 7), (4, 8),
             (5, 9), (6, 10), (7, 10), (7, 11), (8, 12), (9, 13), (10, 14),
             (11, 15), (12, 16), (13, 17), (14, 18), (15, 19), (16, 20),
             (17, 21), (18, 21), (22, 23), (23, 24)]
    for a, b in edges:
        g.add_contact("P%02d" % a, "P%02d" % b)
    return g


def module2(g):
    while True:
        print("\n--- MODULE 2: BFS Contact Tracing ---")
        print("1. Show contact list    2. Trace from index patient (levels)")
        print("3. Shortest chain between two people    4. Add person")
        print("5. Add contact    0. Back")
        ch = input("Choice: ").strip()


        if ch == '0':
            return
        if ch == '1':
            for i in range(len(g.names)):
                print(" %s -> %s" % (g.names[i], ', '.join(g.names[j] for j in g.adj[i]) or '(none)'))
        elif ch == '2':
            s = g.index_of(input("Index patient ID (e.g. P01): ").strip().upper())
            if s == -1:
                print("Unknown person.")
                continue
            level, _ = g.bfs(s)
            mx = max(level)
            for L in range(mx + 1):
                grp = [g.names[i] for i in range(len(level)) if level[i] == L]
                print(" Level %d: %s" % (L, ', '.join(grp)))
            un = [g.names[i] for i in range(len(level)) if level[i] == -1]
            print(" Not connected:", ', '.join(un) if un else 'none')
        elif ch == '3':
            a = input("From: ").strip().upper()
            b = input("To : ").strip().upper()
            if g.index_of(a) == -1 or g.index_of(b) == -1:
                print("Unknown person.")
                continue
            p = g.chain(a, b)
            print(" No chain: not connected." if p is None else
                   " Chain (%d contacts): %s" % (len(p) - 1, ' -> '.join(p)))
        elif ch == '4':
            n = input("New ID: ").strip().upper()
            print(" Added." if n and g.add_person(n) else " Invalid or duplicate ID.")
        elif ch == '5':
            a = input("Person A: ").strip().upper()
            b = input("Person B: ").strip().upper()
            print(" Contact added." if g.add_contact(a, b) else " Invalid, same person, or already exists.")
        else:
            print("Invalid choice.")


# =====================================================================
# MODULE 3 : Recursive vs iterative bacterial growth (Unit 5, T15)
# Model: N(0)=N0, N(t)=r*N(t-1)   (r = growth factor per time step)
# Recurrence for cost: T(n)=T(n-1)+c, T(0)=d -> T(n)=c*n+d = Theta(n)
# (backward substitution: T(n)=T(n-2)+2c = ... = T(0)+n*c)
# =====================================================================
def growth_recursive(n, n0, r, ctr):
    """Counts multiplications and calls. Time Theta(n), space Theta(n) (stack)."""
    ctr[1] += 1                      # call count
    if n == 0:
        return n0
    ctr[0] += 1                      # one multiplication
    return r * growth_recursive(n - 1, n0, r, ctr)


def growth_iterative(n, n0, r, ctr):
    """Time Theta(n), space Theta(1)."""
    val = n0
    for _ in range(n):
        val = val * r
        ctr[0] += 1
    return val


def growth_fast(n, n0, r):
    """Extra: exponentiation by squaring, Theta(log n) multiplications."""
    res, base, e, mults = 1, r, n, 0
    while e > 0:
        if e & 1:
            res *= base
            mults += 1
        base *= base
        mults += 1
        e >>= 1
    return n0 * res, mults


def module3():
    while True:
        print("\n--- MODULE 3: Bacterial Growth (Recursive vs Iterative) ---")
        print("1. Compute colony count      2. Run analysis (4+ sizes, table + plot)")


        print("3. Show recurrence analysis        0. Back")
        ch = input("Choice: ").strip()
        if ch == '0':
              return
        if ch == '1':
              n0 = read_float("Initial colony count N0 (>0): ", 1, 1e9)
              r = read_float("Growth factor per hour r (1 to 10): ", 1, 10)
              n = int(read_float("Hours n (0 to 1000): ", 0, 1000))
              c1, c2 = [0, 0], [0, 0]
              a = growth_recursive(n, n0, r, c1)
              b = growth_iterative(n, n0, r, c2)
              print(" Recursive : %.6g (multiplications=%d, calls=%d)" % (a, c1[0], c1[1]))
              print(" Iterative : %.6g (multiplications=%d)" % (b, c2[0]))
        elif ch == '2':
              sizes = [10, 50, 100, 200, 400, 800]
              print("\n n       rec_mult rec_calls it_mult     rec_time(us) it_time(us)")
              rm, im = [], []
              for n in sizes:
                   c1, c2 = [0, 0], [0, 0]
                   t0 = time.perf_counter()
                   for _ in range(200):
                        c1 = [0, 0]
                        growth_recursive(n, 1, 2, c1)
                   t1 = (time.perf_counter() - t0) / 200 * 1e6
                   t0 = time.perf_counter()
                   for _ in range(200):
                        c2 = [0, 0]
                        growth_iterative(n, 1, 2, c2)
                   t2 = (time.perf_counter() - t0) / 200 * 1e6
                   rm.append(c1[0]); im.append(c2[0])
                   print(" %-5d %-8d %-9d %-9d %-12.2f %.2f" % (n, c1[0], c1[1], c2[0], t1, t2))
              try:
                   import matplotlib
                   matplotlib.use('Agg')
                   import matplotlib.pyplot as plt
                   plt.plot(sizes, rm, 'o-', label='Recursive multiplications')
                   plt.plot(sizes, im, 's--', label='Iterative multiplications')
                   plt.xlabel('Input size n (hours)')
                   plt.ylabel('Operation count')
                   plt.title('Bacterial growth: operations vs n')
                   plt.legend(); plt.grid(True)
                   plt.savefig('growth_operations.png', dpi=150)
                   plt.close()
                   print(" Graph saved as growth_operations.png")
              except ImportError:
                   print(" matplotlib not installed: install it to get the plot (pip install matplotlib)")
        elif ch == '3':
              print("""
 N(t) = r*N(t-1), N(0)=N0
 Cost recurrence : T(n) = T(n-1) + c , T(0) = d
 Backward substitution:
   T(n) = T(n-1)+c = T(n-2)+2c = ... = T(n-k)+k*c ; set k=n
   T(n) = T(0)+n*c = c*n + d => Theta(n)
 Recursive : time Theta(n), space Theta(n) (call stack)
 Iterative : time Theta(n), space Theta(1)
 Best = worst = average case (no input-dependent branches)
 Big-O: O(n)     Omega: Omega(n)      Theta: Theta(n)
 Optional: squaring method needs Theta(log n) multiplications.""")
        else:
              print("Invalid choice.")


# =====================================================================
# TEST SUITE : 3 test cases per module (9 total)
# =====================================================================
def run_tests():
    results = []


    def check(name, expected, actual):
         ok = (expected == actual)
         results.append(ok)
         print("%-4s %-46s expected=%-28s actual=%s" % ("PASS" if ok else "FAIL", name, expected, actual))

    print("\n=== Module 1 tests ===")
    root = build_tree(infix_to_postfix(tokenize("w / h ^ 2")))
    check("T1.1 BMI w=70,h=1.75 (normal)", "22.8571", "%.4f" % evaluate(root, {'w': 70, 'h': 1.75}))


    root = build_tree(infix_to_postfix(tokenize("( s + 2 * d ) / 3")))
    check("T1.2 MAP 120/80 postfix", "s 2 d * + 3 /", ' '.join(infix_to_postfix(tokenize("( s + 2 * d ) / 3"))))
    check("T1.2b MAP 120/80 value", "93.3333", "%.4f" % evaluate(root, {'s': 120, 'd': 80}))
    try:
         r2 = build_tree(infix_to_postfix(tokenize("a / b")))
         evaluate(r2, {'a': 1, 'b': 0})
         got = "no error"
    except ZeroDivisionError:
         got = "ZeroDivisionError"
    check("T1.3 divide by zero (error)", "ZeroDivisionError", got)
    try:
         infix_to_postfix(tokenize("( 1 + 2"))
         got = "no error"
    except ValueError:
         got = "ValueError"
    check("T1.4 mismatched bracket (error)", "ValueError", got)

    print("\n=== Module 2 tests ===")
    g = sample_graph()
    lv, _ = g.bfs(g.index_of("P01"))
    check("T2.1 level of P06 from P01 (normal)", 2, lv[g.index_of("P06")])
    check("T2.2 shortest chain P01->P10", "P01 -> P02 -> P06 -> P10", ' -> '.join(g.chain("P01",
"P10")))
    check("T2.3 P22 unreachable from P01 (boundary)", True, g.chain("P01", "P22") is None)
    check("T2.4 duplicate contact rejected (error)", False, g.add_contact("P01", "P02"))

    print("\n=== Module 3 tests ===")
    c = [0, 0]
    check("T3.1 N0=1,r=2,n=10 recursive", 1024, growth_recursive(10, 1, 2, c))
    c = [0, 0]
    growth_iterative(0, 5, 2, c)
    check("T3.2 n=0 boundary: result & 0 mults", (5, 0), (growth_iterative(0, 5, 2, [0, 0]), c[0]))
    a, b = [0, 0], [0, 0]
    ra = growth_recursive(50, 3, 2, a)
    rb = growth_iterative(50, 3, 2, b)
    check("T3.3 recursive == iterative, n=50", True, ra == rb and a[0] == b[0] == 50)
    check("T3.4 squaring method matches n=30", 1 * 2 ** 30, growth_fast(30, 1, 2)[0])

    print("\nPassed %d of %d" % (sum(results), len(results)))


# =====================================================================
# MAIN MENU
# =====================================================================
def main():
    g = sample_graph()
    while True:
        print("\n========== Group 7 - Biomedical DSA Project ==========")
        print("1. Module 1 - Expression tree (clinical formulas)")
        print("2. Module 2 - BFS contact tracing")
        print("3. Module 3 - Bacterial growth analysis")
        print("4. Run all test cases")
        print("0. Exit")
        ch = input("Choice: ").strip()
        if ch == '1':
            module1()
        elif ch == '2':
            module2(g)
        elif ch == '3':
            module3()
        elif ch == '4':
            run_tests()
        elif ch == '0':
            print("Goodbye.")
            break
        else:
            print("Invalid choice.")


if __name__ == "__main__":
    main()
