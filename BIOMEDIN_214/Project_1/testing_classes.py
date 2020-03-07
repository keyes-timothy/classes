class Shark:
    def __init__(self, name):
        self.name = name
        
    def swim(self):
        # Reference the name
        print(self.name + " is swimming.")
        
    def be_awesome(self):
        # Reference the name
        print(self.name + " is being awesome.")


sammy = Shark("Sammy")
sammy.swim()
sammy.be_awesome()

