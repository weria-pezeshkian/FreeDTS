/*
 ===============================================================
 Generic Factory (Templated)
 ===============================================================
 Weria Pezeshkian (2026)
 
 This class implements a reusable, type-safe singleton factory
 for creating objects derived from a specified base class `Base`.

 The factory enables dynamic runtime construction of derived
 objects using a string identifier, without requiring direct
 knowledge of the concrete types at the call site.

 ---------------------------------------------------------------
 Key Design Idea
 ---------------------------------------------------------------
 The factory is templated on the base class:

        Factory<Base>

 This allows a single implementation to be reused for multiple
 independent class hierarchies, such as:

        using FactoryDynamicBox =
            Factory<AbstractDynamicBox>;


 Each specialization maintains its own registry and singleton
 instance.

 ---------------------------------------------------------------
 Creator Function Interface
 ---------------------------------------------------------------
 Each derived class must provide (or be wrapped by) a creator
 function with the following signature:

        Base* Create(std::istream& input, State* state);

 The creator is responsible for:
   - Parsing its own parameters from the input stream
   - Constructing the derived object
   - Passing arguments to the constructor in the correct order

 NOTE:
 The factory does NOT enforce a specific constructor signature,
 but the creator must adapt whatever constructor is used.

 ---------------------------------------------------------------
 How It Works
 ---------------------------------------------------------------
 1. A derived class defines a creator function (or lambda):

        static Base* Create(std::istream& input, State* state);

 2. The creator is registered with a string key:

        Factory<Base>::Instance().Register("TypeName", Create);

 3. A caller requests an object:

        Base* obj =
            Factory<Base>::Instance()
                .Create(name, input, state);

 4. The factory:
      - Looks up the string identifier
      - Calls the associated creator
      - Returns the constructed object

 ---------------------------------------------------------------
 Parameters
 ---------------------------------------------------------------
 name  : string identifier of the desired derived type
 input : stream containing parameters for construction
 state : pointer to State object

 The `state` pointer is always forwarded to the creator, but its
 position in the derived constructor is not dictated by the
 factory.

 ---------------------------------------------------------------
 Return Value
 ---------------------------------------------------------------
 Returns:
    Pointer to a newly created object of a type derived from Base

 Returns nullptr if:
    - The name is not registered
    - The creator fails or returns nullptr

 ---------------------------------------------------------------
 Lifetime & Ownership
 ---------------------------------------------------------------
 The factory does NOT manage memory.

 The caller is responsible for:
    - Owning the returned pointer
    - Deleting it when appropriate


 ---------------------------------------------------------------
 Registration Notes
 ---------------------------------------------------------------
 - Registration is typically done at static initialization time
   inside a .cpp file.

 - Multiple factories (different Base types) are completely
   independent due to template instantiation.

 ---------------------------------------------------------------
 Thread Safety
 ---------------------------------------------------------------
 - Singleton initialization is thread-safe
 - Registration and creation are NOT inherently thread-safe
   and must be externally synchronized if used concurrently

 ---------------------------------------------------------------
 Summary
 ---------------------------------------------------------------
 This templated factory removes duplication across multiple
 factory implementations while preserving flexibility in how
 derived objects are constructed.

 It provides:
   - Runtime polymorphic construction
   - Decoupling between caller and concrete types
   - Reusable infrastructure across class hierarchies

 ===============================================================
*/
template <typename Base>
class Factory {
public:
    using Creator = std::function<Base*(std::istream&, State*)>;

    static Factory& Instance()
    {
        static Factory instance;
        return instance;
    }

    void Register(const std::string& name, Creator c)
    {
        m_Creators[name] = std::move(c);
    }

    Base* Create(
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
    Factory() = default;
    Factory(const Factory&) = delete;
    Factory& operator=(const Factory&) = delete;

    std::unordered_map<std::string, Creator> m_Creators;
};


// second type of factory for global constraint
class VAHGlobalMeshProperties;
template <typename Base>
class FactoryGlobal {
public:
    using Creator = std::function<Base*(std::istream&, VAHGlobalMeshProperties*)>;

    static FactoryGlobal& Instance()
    {
        static FactoryGlobal instance;
        return instance;
    }

    void Register(const std::string& name, Creator c)
    {
        m_Creators[name] = std::move(c);
    }

    Base* Create(
        const std::string& name,
        std::istream& input,
        VAHGlobalMeshProperties* VAH)
    {
        auto it = m_Creators.find(name);

        if (it == m_Creators.end())
            return nullptr;

        return it->second(input, VAH);
    }

private:
    FactoryGlobal() = default;
    FactoryGlobal(const FactoryGlobal&) = delete;
    FactoryGlobal& operator=(const FactoryGlobal&) = delete;

    std::unordered_map<std::string, Creator> m_Creators;
};
