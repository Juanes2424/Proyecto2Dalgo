from collections import deque, defaultdict
import random
import heapq
import time


def bfs_capacity_path(graph, source, sink, parent, excluded_node):
    """Breadth-first search to find an augmenting path in the residual graph, excluding a node."""
    visited = set()
    queue = deque([source])
    visited.add(source)

    while queue:
        u = queue.popleft()

        for v in graph[u]:
            if v == excluded_node:  # Skip the excluded node
                continue
            if (
                v not in visited and graph[u][v] > 0
            ):  # Check if not visited and capacity exists
                parent[v] = u
                visited.add(v)
                if v == sink:
                    return True
                queue.append(v)

    return False


def EK(graph, source, sink, excluded_node=None):
    """Implementation of the Edmonds-Karp algorithm with the option to exclude a node.
    Returns the flow for each edge along with the maximum flow."""
    max_flow = 0
    parent = {}
    flow_graph = {u: {v: 0 for v in graph[u]} for u in graph}  # Initialize flow graph

    # Run BFS as long as there's a path from source to sink
    while bfs_capacity_path(graph, source, sink, parent, excluded_node):
        # Find the maximum flow through the path found by BFS
        path_flow = float("Inf")
        s = sink
        while s != source:
            path_flow = min(path_flow, graph[parent[s]][s])
            s = parent[s]

        # Update residual capacities of the edges and reverse edges
        v = sink
        while v != source:
            u = parent[v]
            # Update the capacity of the forward edge
            graph[u][v] -= path_flow
            flow_graph[u][v] += path_flow  # Record flow on the edge

            # Ensure the reverse edge exists and update its capacity
            if v not in graph:
                graph[v] = {}
            if u in graph[v]:
                graph[v][u] += path_flow
            else:
                # Initialize the reverse edge with path_flow if it doesn't exist
                graph[v][u] = path_flow

            v = parent[v]

        max_flow += path_flow

    return flow_graph, max_flow


def bruteForce():
    _, maxFlow = EK(shallow_copy_graph(grafo), -1, -2, None)

    nodo = None
    minSoFar = maxFlow

    for i in tipos[1]:
        _, res = EK(shallow_copy_graph(grafo), -1, -2, i)
        if res < minSoFar:
            minSoFar = res
            nodo = i

    # print(nodo, maxFlow, minSoFar)

    return nodo, maxFlow, minSoFar


def shallow_copy_graph(graph):
    # This function will make a shallow copy of the graph, including deep copies of nested dictionaries
    return {k: v.copy() if isinstance(v, dict) else v for k, v in graph.items()}


def greedy():
    nodo = None
    maxSoFar = -1

    for i in tipos[1]:
        if min(entrantes[i], salientes[i]) >= maxSoFar:
            maxSoFar = min(entrantes[i], salientes[i])
            nodo = i

    return (
        nodo,
        EK(shallow_copy_graph(grafo), -1, -2, None),
        EK(shallow_copy_graph(grafo), -1, -2, nodo),
    )


def greedy2():
    g, mF = EK(shallow_copy_graph(grafo), -1, -2, None)

    nodoMax = None
    maxFlowPerNode = -1

    for llave in g:
        if llave not in tipos[1]:
            continue

        cou = 0
        for v, peso in g[llave].items():
            cou += peso

        if cou > maxFlowPerNode:
            maxFlowPerNode = cou
            nodoMax = llave

    # print(nodoMax, maxFlowPerNode)
    # print(g)

    chosenNode, minFlow = None, mF

    for v in g[nodoMax]:
        if v in tipos[1]:
            # print("V " + str(v))
            _, flujo = EK(shallow_copy_graph(grafo), -1, -2, v)

            if flujo < minFlow:
                chosenNode = v
                minFlow = flujo

    # print("V " + str(nodoMax))
    _, flujo = EK(shallow_copy_graph(grafo), -1, -2, nodoMax)

    if flujo < minFlow:
        chosenNode = nodoMax
        minFlow = flujo

    # print(chosenNode, mF, minFlow)

    return chosenNode, mF, minFlow


def greedy3():

    g, maxFlow = EK(shallow_copy_graph(grafo), -1, -2, None)

    # print(g)

    pq = []

    print("IG")

    for i in tipos[1]:
        cou = 0
        for j, w in g[i].items():
            cou += w

        pq.append((-cou, i))

    print("HOLA")

    pq.sort()

    elemento_1 = pq[0]

    nodoAQuitar = elemento_1[1]
    minFlow = EK(shallow_copy_graph(grafo), -1, -2, nodoAQuitar)[1]

    print("EMPATE")

    indice = 0

    start_time = time.time()

    while pq:
        print(indice)

        if time.time() - start_time > 10:
            print("Time limit exceeded, moving on...")
            break

        ele = pq[indice]
        if ele[0] == elemento_1[0]:
            nuevoMin = EK(shallow_copy_graph(grafo), -1, -2, ele[1])[1]
            if nuevoMin < minFlow:
                nodoAQuitar = ele[1]
                minFlow = nuevoMin
        else:
            break  # Stop if priorities no longer match

        indice += 1

    return (nodoAQuitar, maxFlow, minFlow)


def greedy4():
    pass


def procesarEntrada():
    primera_linea = input().split(" ")
    n, d = int(primera_linea[0]), int(primera_linea[1])

    for _ in range(n):
        celula = input().strip().split()
        id, x, y, tipo, peptidos = (
            int(celula[0]),
            int(celula[1]),
            int(celula[2]),
            int(celula[3]),
            celula[4:],
        )

        peptidos.sort()

        grafo[id] = {}
        entrantes[id] = 0
        salientes[id] = 0

        entradas[id] = (x, y, peptidos)
        tipos[tipo - 1].append(id)

    return d


def crearGrafo():
    # conexiones calculadoras
    for i in range(len(tipos[1])):
        for j in range(i + 1, len(tipos[1])):
            u, v = tipos[1][i], tipos[1][j]
            w = numeroMensajes(u, v)

            if w != 0:
                entrantes[u] = entrantes[u] + w
                salientes[u] = salientes[u] + w
                entrantes[v] = entrantes[v] + w
                salientes[v] = salientes[v] + w

                grafo[u].update({v: w})
                grafo[v].update({u: w})

    for i in tipos[0]:
        for j in tipos[1]:
            w = numeroMensajes(i, j)

            if w != 0:
                salientes[i] = salientes[i] + w
                entrantes[j] = entrantes[j] + w
                grafo[i].update({j: w})

    for i in tipos[1]:
        for j in tipos[2]:
            w = numeroMensajes(i, j)

            if w != 0:
                salientes[i] = salientes[i] + w
                entrantes[j] = entrantes[j] + w
                grafo[i].update({j: w})


def numeroMensajes(m, n):
    a, b = entradas[m], entradas[n]

    if (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 > d**2:
        return 0

    i, j = 0, 0
    cou = 0

    while i < len(a[2]) and j < len(b[2]):
        if a[2][i] == b[2][j]:
            cou += 1
            i += 1
            j += 1
        elif a[2][i] < b[2][j]:
            i += 1
        else:
            j += 1

    return cou


def crearSinkSource():
    grafo[-1] = {}  # sourc
    grafo[-2] = {}  # sink

    for i in tipos[0]:
        grafo[-1].update({i: 10**8})

    for i in tipos[2]:
        grafo[i].update({-2: 10**8})


def simularEntrada():
    n_nodos, d = random.randint(1000, 1001), random.randint(1, 3)
    indice = 1

    for i in range(0, int(n_nodos ** (1 / 2) + 2)):
        for j in range(0, int(n_nodos ** (1 / 2) + 2)):
            grafo[indice] = {}
            entrantes[indice] = 0
            salientes[indice] = 0

            entradas[indice] = (i, j, ["abcde"])
            tipos[(i) % 3].append(indice)

            indice += 1

    # print(grafo)

    return d


"""
for _ in range(10):
    tipos = [[], [], []]
    entradas = {}
    grafo = {}

    entrantes = {}
    salientes = {}

    d = simularEntrada()
    crearGrafo()
    crearSinkSource()

    try:

        a = greedy3()

        # print("GREEDY")
        # greedy()

        # print(a)
        print("AHM")

        b = bruteForce()

        if a != b:
            print(">>>>>>>>>>>>>>WRONG ")
        else:
            print(">>>>>>>>>>>>>>>>>ALL GOOD")
    except:
        print("0 0 0")

"""

casos = int(input())

for _ in range(casos):

    try:
        tipos = [[], [], []]
        entradas = {}
        grafo = {}

        entrantes = {}
        salientes = {}

        d = procesarEntrada()
        crearGrafo()
        crearSinkSource()

        # print(grafo)
        print(grafo)

        print("GREEDY 3")

        print()
        print()
        print()

        print(greedy3())
        print(bruteForce())

        print()
        print()
        print()

        # print("GREEDY")
        # greedy()

        # print("BRUTE FORCE")
        # print(bruteForce())
    except:
        print("0 0 0")
