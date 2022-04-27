import math
import random
from matplotlib import pyplot as plt
import pandas as pd

class Edge:
    def __init__(self, a, b, edge_distance, initial_pheromone):
        self.a = a
        self.b = b
        self.edge_distance = edge_distance
        self.pheromone = initial_pheromone

class Ant:
    def __init__(self, alpha, beta, num_nodes, dist_matrix):
        self.alpha = alpha
        self.beta = beta
        self.num_nodes = num_nodes
        self.dist_matrix = dist_matrix
        self.tour = None
        self.distance = 0.0

    def probability_to_node(self,unvisited_node):
      dist=pow(self.dist_matrix[self.tour[-1]][unvisited_node].pheromone,\
               self.alpha)*pow((1/self.dist_matrix[self.tour[-1]]\
                                [unvisited_node].edge_distance),self.beta)
      return dist

    def roulette_wheel(self,unvisited_nodes):
        roulette_wheel = 0.0
        for unvisited_node in unvisited_nodes:
            roulette_wheel += self.probability_to_node(unvisited_node)

        random_value = random.uniform(0.0, roulette_wheel)
        wheel_position = 0.0

        for unvisited_node in unvisited_nodes:
            wheel_position += self.probability_to_node(unvisited_node)
            if wheel_position >= random_value:
                return unvisited_node
            
    def calculate_tour_distance(self):
        self.distance = 0.0
        for i in range(len(self.tour)-1):
            self.distance +=\
            self.dist_matrix[self.tour[i]][self.tour[(i + 1)]].edge_distance
        return self.distance

class ACO:
    def __init__(self, mode='ACO', colony_size=10, alpha=1.5, beta=3.0,
                 rho=0.1, initial_pheromone=1.0, steps=100, nodes=None,\
                 labels=None):
        self.mode = mode
        self.colony_size = colony_size
        self.rho = rho           
        self.steps = steps
        self.num_nodes = len(nodes)
        self.nodes = nodes
        self.global_best_tour = None
        self.global_best_distance = float("inf")
        self.dist_matrix = self.get_distance_matrix(initial_pheromone)
        if labels is not None:
            self.labels = labels
        else:
            self.labels = range(0, self.num_nodes + 0)

    def get_distance_matrix(self,initial_pheromone):
        dist_matrix = [[None] * self.num_nodes for i in range(self.num_nodes)]
        for i in range(self.num_nodes):
            for j in range(i + 1, self.num_nodes):
              dist_matrix[i][j]=self.Edge(i, j,self.get_distance_two_nodes(i, j),
                                                            initial_pheromone)
              dist_matrix[j][i]=dist_matrix[i][j]
        return dist_matrix

    def get_distance_two_nodes(self,i,j):
        distance = math.sqrt(pow(self.nodes[i][0] - self.nodes[j][0], 2.0) + 
                             pow(self.nodes[i][1] - self.nodes[j][1], 2.0))        
        return distance

    def add_pheromone(self, tour, distance, pheromone_factor=1.0):
        pheromone_to_add = self.rho/ (self.num_nodes / distance)
        self.local_pheromone_update(tour,pheromone_to_add)

    def local_pheromone_update(self,tour,pheromone_to_add):
        for i in range(len(tour)-1):
            self.dist_matrix[tour[i]][tour[(i + 1)]].pheromone +=pheromone_to_add

    def pheromone_disintegration(self):
        for i in range(self.num_nodes):
            for j in range(i + 1, self.num_nodes):
                self.dist_matrix[i][j].pheromone *= (1.0 - self.rho)

    def global_pheromone_update(self):
        for i in range(len(self.global_best_tour)-1):
         self.dist_matrix[self.global_best_tour[i]]\
                         [self.global_best_tour[(i+1)]].pheromone\
                         += self.rho*(1/self.global_best_distance)

    def aco(self):
        for step in range(self.steps):
            for ant in self.ants:
                self.add_pheromone(ant.find_tour(), ant.calculate_tour_distance())
                if ant.distance < self.global_best_distance:
                    self.global_best_tour = ant.tour
                    self.global_best_distance = ant.distance
                self.pheromone_disintegration()

                self.global_pheromone_update()

    def run(self):
     print('Started')
     if self.mode == 'ACO':
       self.aco()
       print('Ended')
       print('Sequence:<-{0}->'.format('-'.join(str(self.labels[i])\
             for i in self.global_best_tour)))
       print('Total distance travelled:{0}\n'.format\
             (round(self.global_best_distance,2)))

    def plot(self, line_width=1, point_radius=math.sqrt(2.0),
             annotation_size=8, dpi=120, save=True, name=None):
        x = [self.nodes[i][0] for i in self.global_best_tour]
        y = [self.nodes[i][1] for i in self.global_best_tour]
        plt.plot(x, y, linewidth=line_width)
        plt.scatter(x, y, s=math.pi * (point_radius ** 2.0))
        plt.title(self.mode)
        for i in self.global_best_tour:
            plt.annotate(self.labels[i], self.nodes[i], size=annotation_size)
        if save:
            if name is None:
                name = '{0}.png'.format(self.mode)
            plt.savefig(name, dpi=dpi)
        plt.show()
        plt.gcf().clear()

class TSP(ACO):      #Traveling salesman problem
    class Edge(Edge):
        def __init__(self, a, b, edge_distance, initial_pheromone):
            super().__init__(a, b, edge_distance, initial_pheromone)
    class Ant(Ant):
        def __init__(self, alpha, beta, num_nodes, dist_matrix):
            super().__init__(alpha, beta, num_nodes, dist_matrix)

        def select_node(self):
            unvisited_nodes = [node for node in range(self.num_nodes)\
                               if node not in self.tour]
            selected_node=self.roulette_wheel(unvisited_nodes)
            return selected_node

        def find_tour(self):
            self.tour = [random.randint(0, self.num_nodes - 1)]
            while len(self.tour) < self.num_nodes:
                self.tour.append(self.select_node())
            self.tour.append(self.tour[0])
            return self.tour

    def __init__(self,mode='ACO',colony_size=10,alpha=1.5,beta=2.0,
                 rho=0.1,initial_pheromone=100.0,steps=100,nodes=None,\
                 labels=None):
        
        super().__init__(mode, colony_size, alpha, beta,
                      rho,initial_pheromone, steps, nodes, labels)

        self.ants = [self.Ant(alpha, beta, self.num_nodes, self.dist_matrix)\
                    for _ in range(self.colony_size)]


class CVRP(ACO):     #Capacitated vehicle routing problem VRP
    class Node:
        def __init__(self,weight):
            self.weight=weight
    class Edge(Edge):
        def __init__(self, a, b, edge_distance, initial_pheromone):
            super().__init__(a, b, edge_distance, initial_pheromone)
    class Ant(Ant):
        def __init__(self, alpha, beta, num_nodes, dist_matrix,
                     weight_matrix,vehicle_num,load_capacity):
            super().__init__(alpha, beta, num_nodes, dist_matrix)
            self.weight_matrix=weight_matrix
            self.current_load=0.0
            self.load_capacity=load_capacity
            self.vehicle_num=vehicle_num
            self.current_vehicle_num=vehicle_num
        def select_node(self):
            unvisited_nodes = [node for node in range(self.num_nodes)\
                               if node not in self.tour]
            unvisited_nodes_selection=[ ]
            for next_node in unvisited_nodes:
                if(self.current_load + self.weight_matrix[next_node].weight>\
                   self.load_capacity):
                    continue
                unvisited_nodes_selection.append(next_node)
            if not unvisited_nodes_selection:
                self.current_load=0;
                self.current_vehicle_num-=1
                return 0;
            selected_node=self.roulette_wheel(unvisited_nodes_selection)
            self.current_load+=self.weight_matrix[selected_node].weight
            return selected_node

        def find_tour(self):
            self.current_vehicle_num=self.vehicle_num
            self.current_load=0
            self.tour =[0]
            while len(self.tour) < self.num_nodes+self.vehicle_num :
                self.tour.append(self.select_node())
            if(self.current_vehicle_num < 0 or self.tour[-1]!=0 ):
                return self.find_tour() 
            return self.tour

    def __init__(self, mode='ACO', colony_size=10, alpha=1.5, beta=2.5,
                 rho=0.1,initial_pheromone=100, steps=100, nodes=None,
                 nodes_weight=None,vehicle_num=4,vehicle_capacity=15,labels=None):
      super().__init__(mode, colony_size, alpha, beta,
                      rho,initial_pheromone,steps,nodes,labels)
      self.nodes_weight=self.get_weight_matrix(nodes_weight)
      self.vehicle_num=vehicle_num
      self.vehicle_capacity=vehicle_capacity
      self.ants=[self.Ant(alpha, beta, self.num_nodes, self.dist_matrix,\
                        self.nodes_weight,self.vehicle_num,\
                        self.vehicle_capacity)for _ in range(self.colony_size)]
 
    def get_weight_matrix(self,nodes_weight):
        weight_matrix=[[None] * self.num_nodes for i in range(self.num_nodes)]
        for i in range (self.num_nodes):
            weight_matrix[i]=self.Node(self.get_edge_weight(nodes_weight,i))
        return weight_matrix

    def get_edge_weight(self,nodes_weight,i):
        weight=nodes_weight[i]
        return weight

if __name__ == '__main__':

    _colony_size = 100
    _steps = 30


    data=pd.read_csv("matrix_weight_4_vehicles.csv")
    lat=data['lat'].tolist()
    long=data['long'].tolist()
    weight_from_csv=data['weight'].tolist()
    array_from_csv=list(zip(lat,long))

    acs = TSP(mode='ACO', colony_size=_colony_size, steps=_steps, nodes=array_from_csv)
    vrp = CVRP(mode='ACO', colony_size=_colony_size, steps=_steps, nodes=array_from_csv,nodes_weight=weight_from_csv)


    acs.run()
    acs.plot()

    vrp.run()
    vrp.plot()





































