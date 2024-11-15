# Autores:
# Juan Esteban Rios Gonzalez 202221404
# Jessica Sofia Garay Acosta 202310514

from collections import deque
import time


def bfs(graph, source, sink, parent, nodoExcluido):
    visited = set()
    q = deque([source])
    visited.add(source)

    while q:
        u = q.popleft()

        for v in graph[u]:
            if v == nodoExcluido:
                continue

            if v not in visited and graph[u][v] > 0:
                parent[v] = u
                visited.add(v)

                if v == sink:
                    return True

                q.append(v)

    return False


def EK(graph, source, sink, excluded_node=None):
    maxFlow = 0
    parent = {}
    grafoFlujos = {u: {v: 0 for v in graph[u]} for u in graph}

    while bfs(graph, source, sink, parent, excluded_node):
        flujoDeCamino = float("Inf")
        s = sink

        while s != source:
            flujoDeCamino = min(flujoDeCamino, graph[parent[s]][s])
            s = parent[s]

        v = sink

        while v != source:
            u = parent[v]
            graph[u][v] -= flujoDeCamino
            grafoFlujos[u][v] += flujoDeCamino

            if v not in graph:
                graph[v] = {}

            if u in graph[v]:
                graph[v][u] += flujoDeCamino
            else:
                graph[v][u] = flujoDeCamino

            v = parent[v]

        maxFlow += flujoDeCamino

    return grafoFlujos, maxFlow


def copiaGrafo(graph):
    return {k: v.copy() if isinstance(v, dict) else v for k, v in graph.items()}


def greedy3():

    g, maxFlow = EK(copiaGrafo(grafo), -1, -2, None)
    pq = []

    for i in tipos[1]:
        cou = 0
        for j, w in g[i].items():
            cou += w

        pq.append((-cou, i))

    pq.sort()
    elemento_1 = pq[0]

    nodoAQuitar = elemento_1[1]
    minFlow = EK(copiaGrafo(grafo), -1, -2, nodoAQuitar)[1]

    indice = 0

    start_time = time.time()

    while pq:
        if time.time() - start_time > 10:
            break

        if indice >= len(pq):
            break

        ele = pq[indice]
        if ele[0] == elemento_1[0]:
            nuevoMin = EK(copiaGrafo(grafo), -1, -2, ele[1])[1]

            if nuevoMin < minFlow:
                nodoAQuitar = ele[1]
                minFlow = nuevoMin
        else:
            break

        indice += 1

    return (nodoAQuitar, maxFlow, minFlow)


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
                grafo[u].update({v: w})
                grafo[v].update({u: w})

    for i in tipos[0]:
        for j in tipos[1]:
            w = numeroMensajes(i, j)

            if w != 0:
                grafo[i].update({j: w})

    for i in tipos[1]:
        for j in tipos[2]:
            w = numeroMensajes(i, j)

            if w != 0:
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
    grafo[-1] = {}  # source
    grafo[-2] = {}  # sink

    for i in tipos[0]:
        grafo[-1].update({i: 10**8})

    for i in tipos[2]:
        grafo[i].update({-2: 10**8})


casos = int(input())

for _ in range(casos):

    try:
        tipos = [[], [], []]
        entradas = {}
        grafo = {}

        d = procesarEntrada()
        crearGrafo()
        crearSinkSource()

        respuesta = greedy3()

        print(str(respuesta[0]) + " " + str(respuesta[1]) + " " + str(respuesta[2]))
    except:
        print("0 0 0")
