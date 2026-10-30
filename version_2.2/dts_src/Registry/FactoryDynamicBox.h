/*
 ===============================================================
 FactoryDynamicBox
 ===============================================================

 This class implements a singleton factory for creating instances of
 `AbstractVertexPositionIntegrator` derived classes dynamically at runtime
 based on a string identifier.

 IMPORTANT:
 The factory does NOT enforce any specific constructor signature
 for derived classes.

 Each derived class is responsible for:
   - Parsing its own parameters from the input stream
   - Calling its constructor with the correct argument order

 The only required interface is the creator function.

 ---------------------------------------------------------------
 How it works:
 ---------------------------------------------------------------
 1. Each derived class (or its .cpp file) defines a creator function:

        static AbstractVertexPositionIntegrator* Create(
            std::istream& input,
            State* state);

    or an equivalent free/static function.

 2. The creator:
      - Reads parameters from the input stream
      - Constructs the derived object
      - Passes the State pointer in whatever position the constructor expects

 3. The creator is registered with the factory using a string identifier.

 4. The State class calls:

        FactoryDynamicBox::Instance()
            .Create(name, input, this);

 5. The factory:
      - Looks up the registered creator
      - Calls it
      - Returns the constructed object (or nullptr if not found)

 ---------------------------------------------------------------
 Notes
 ---------------------------------------------------------------
 - The State pointer is provided to the creator, but its position
   in the derived class constructor is NOT fixed by the factory.

 - Returning nullptr indicates:
      * Unknown type name
      * Invalid or failed parameter parsing

 - Registration typically happens at static initialization time
   inside a .cpp file.

 - The factory owns NO memory. The caller is responsible for
   managing (and deleting) the returned object.
 ---------------------------------------------------------------
*/
class State;   // forward declaration
class FactoryDynamicBox {
public:

    /*
    ---------------------------------------------------------------
    Creator Function Type
    ---------------------------------------------------------------
    Every creator must accept:

        std::istream& input
        State* state

    The state pointer is always passed LAST so derived constructors
    follow the rule:

        Derived(...parameters..., State* state)
    */
    using Creator =
        std::function<AbstractDynamicBox*(std::istream&, State*)>;

    /*
    ---------------------------------------------------------------
    Singleton Instance
    ---------------------------------------------------------------
    */
    static FactoryDynamicBox& Instance()
    {
        static FactoryDynamicBox instance;
        return instance;
    }

    /*
    ---------------------------------------------------------------
    Register
    ---------------------------------------------------------------
    Registers a creator function with a string identifier.
    */
    void Register(const std::string& name, Creator c)
    {
        m_Creators[name] = c;
    }

    /*
    ---------------------------------------------------------------
    Create
    ---------------------------------------------------------------
    Creates a derived evolution method.

    Parameters:
        name  : identifier of derived class
        input : stream containing parameters
        state : pointer to the State object creating the evolution

    Returns:
        pointer to AbstractVertexPositionIntegrator or nullptr
    */
    AbstractDynamicBox* Create(
        const std::string& name,
        std::istream& input,
        State* state)
    {
        auto it = m_Creators.find(name);

        if (it == m_Creators.end())
            return nullptr;

        return it->second(input, state);
    }

private:

    FactoryDynamicBox() = default;

    FactoryDynamicBox(
        const FactoryDynamicBox&) = delete;

    FactoryDynamicBox& operator=(
        const FactoryDynamicBox&) = delete;

    std::unordered_map<std::string, Creator> m_Creators;
};
