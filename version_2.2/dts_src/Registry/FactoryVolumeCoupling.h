/*
===============================================================
FactoryVolumeCoupling
===============================================================

This class implements a **singleton factory** for creating instances of
`AbstractVolumeCoupling` derived classes (e.g., Apply_Osmotic_Pressure,
VolumeCouplingSecondOrder, etc.) dynamically at runtime based on a string
identifier.

The factory allows decoupling of object creation from the rest of the code,
enabling polymorphic instantiation without hardcoding derived types.
This is particularly useful when reading simulation configurations from files
and creating the appropriate volume coupling objects on-the-fly.

---------------------------------------------------------------
How it works:
---------------------------------------------------------------
1. Each derived class defines a static "creator" function that reads
   the required parameters from an input stream and constructs an instance.
2. Each derived class registers itself with the factory using its default
   name and creator function. Typically, this is done using a static object
   at global scope.
3. When the simulation code needs a volume coupling object, it calls
   FactoryVolumeCoupling::Instance().Create(name, input, vah),
   which:
   - Looks up the name in the registry.
   - If found, invokes the creator function to produce a new object.
   - Returns a pointer to the newly created AbstractVolumeCoupling.
---------------------------------------------------------------
*/

class FactoryVolumeCoupling {
public:
    // Define the type of the creator function
    // It takes:
    //   - std::istream& input: input stream to read parameters
    //   - VAHGlobalMeshProperties* vah: pointer to shared mesh properties
    // Returns a pointer to a newly constructed AbstractVolumeCoupling
    using Creator = std::function<AbstractVolumeCoupling*(std::istream&, VAHGlobalMeshProperties*)>;

    /*
    ---------------------------------------------------------------
    Instance
    ---------------------------------------------------------------
    Returns a reference to the singleton instance of the factory.
    The static local variable ensures:
      - Only one instance exists during program execution.
      - Thread-safe initialization in C++11 and later.
    */
    static FactoryVolumeCoupling& Instance() {
        static FactoryVolumeCoupling instance;
        return instance;
    }

    /*
    ---------------------------------------------------------------
    Register
    ---------------------------------------------------------------
    Registers a new creator function with a given string identifier.

    Parameters:
      - name: The string name used to identify this derived class.
      - c: The creator function responsible for constructing objects of
           this type.

    Example:
      FactoryVolumeCoupling::Instance().Register(
          "OsmoticPressure",
          ApplyOsmoticPressureRegister::Create
      );
    */
    void Register(const std::string& name, Creator c) {
        m_Creators[name] = c;
    }

    /*
    ---------------------------------------------------------------
    Create
    ---------------------------------------------------------------
    Creates a new AbstractVolumeCoupling object based on a string name.

    Parameters:
      - name: The string identifier for the desired derived class.
      - input: Input stream from which to read constructor parameters.
      - vah: Pointer to shared VAHGlobalMeshProperties.

    Returns:
      - A pointer to a newly created AbstractVolumeCoupling, or nullptr
        if the name is not registered.

    Notes:
      - The caller is responsible for managing the lifetime of the returned
        object (typically using smart pointers is recommended).
      - This function decouples the code that creates objects from the
        knowledge of derived classes, enabling dynamic and extensible
        design.
    */
    AbstractVolumeCoupling* Create(
        const std::string& name,
        std::istream& input,
        VAHGlobalMeshProperties* vah)
    {
        auto it = m_Creators.find(name);
        if(it == m_Creators.end())
            return nullptr; // name not found in registry

        // Call the registered creator function to construct the object
        return it->second(input, vah);
    }

private:
    // Private constructor to enforce singleton pattern
    FactoryVolumeCoupling() = default;

    // Disable copy constructor and assignment operator
    FactoryVolumeCoupling(const FactoryVolumeCoupling&) = delete;
    FactoryVolumeCoupling& operator=(const FactoryVolumeCoupling&) = delete;

    // Registry mapping string names to creator functions
    std::unordered_map<std::string, Creator> m_Creators;
};
