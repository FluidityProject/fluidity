/*
 * optparse.h
 *
 * An option parser for C++ modeled after the optparse Python
 * module. Create an instance of OptionParser in your program
 * and you can add options to it then call the method
 * parse_args with argc and argv and it will parse them. The
 * '-h' and '--help' options are built in, they print the
 * usage message.
 *
 * W. Evan Sheehan <evan_sheehan@nrel.gov>
 *
 */

#ifndef OPTPARSE_H
#define OPTPARSE_H

#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>

namespace optparse {

/**
 * Action type. 
 * Use this type to determine what we're storing, naming
 * convention is the same as the optparse python module.
 * Even with STORE_FALSE and STORE_TRUE the values are
 * going to be stored as strings, false is "0" true is "1"
 * so you can just atoi() the value and use that for tests.
 */
typedef enum {STORE_FALSE=0, STORE_TRUE, STORE} action_t;

/** 
 * Option argument type.
 * Defines expected option argument type.
 */
typedef enum {STRING=0, BOOL, INT, DOUBLE} type_t;

/** Option for OptionParser */
class Option {
public:
    Option (std::string shrt, std::string lng, std::string dest,
            std::string hlp, action_t act, std::string dfault, type_t type,
            std::string allowed_args);
    ~Option ();
    
    /**
     * Validate option argument.
     * Test if argument is valid for the option. The argument must be valid w.r.t. 
     * to type_ and allowed_values_.
     * @param argument Argument to be tested.
     * @return true if argument is valid.
     */
    bool is_allowed(std::string argument);

    action_t action;            // define what we store
    std::string shrt_flag;      // short flag, something like -o
    std::string lng_flag;       // long flag, something like --outfile
    std::string help;           // help message about the option
    std::string destination;    // the key used to store this option in the parser.
    /**
     * Default value. An empty string means no default value exist.
     */
    std::string dfault_;
    /**
     * Expected type the argument of this option should have. 
     * Only relevant for STORE action.
     * Used for generating the usage string.
     */
    type_t type_;
    std::string allowed_args_;        //!< Comma-separated list of allowed option arguments.
};

class OptionError : public std::runtime_error {
public:
    OptionError(std::string option, std::string error_msg)
            : std::runtime_error("") {
        std::ostringstream str;
        str << "OptionParser error: option: " << option << ", cause: " << error_msg;
        msg_ = str.str();
    }
    virtual ~OptionError() throw() {}
    /**
     * Human-readable error message.
     * @return Pointer to error message.
     */
    virtual const char* what() const throw() { return msg_.c_str(); }
protected:
    std::string msg_;
};

/** Option parser. */
class OptionParser {
public:
    OptionParser (std::string usage="");
    virtual ~OptionParser ();

    /**
     * Add an option to the parser.
     * @param shrt_flag Short option name, like e.g. "-q".
     * @param lng_flag Long option name, like e.g. "--quiet"
     * @param destination Key under which the option argument is stored in the dictionary.
     * @param help Help string for generating the usage info.
     * @param act Action, one of STORE, STORE_TRUE, STORE_FALSE.
     * @param type Type info of the expected option argument. One of INT, DOUBLE, STRING, BOOL.
     * @param dfault Default value. Value stored in the dictionary if the option is not given.
     * @param allowed_values List of possible option values. 
     *                       A string of comma-separated allowed values. An empty string means
     *                       that any value is allowed.
     */
    void add_option (std::string shrt_flag, std::string lng_flag,
                     std::string destination, std::string help="", action_t act=STORE,
                     type_t type=STRING, std::string dfault="",
                     std::string allowed_values="");

    /**
     * Get option argument.
     * @param option_name Name of option, as given in "destination" parameter of "add_option".
     * @return Option argument with type appendix removed.
     */
    std::string get_option(std::string option_name);
    template<typename T> T get(std::string option_name);
    /* Parse the commandline args. */
    void parse_args (int argc, char **argv);
    /**
     * Write usage info to stream.
     * The usage info includes a formatted list of all options the parser knows about, 
     * including the help string, expected argument type and default value.
     * @param os Output stream the usage info is written to.
     */
    void help (std::ostream& os);
    /**
     * Convert type_t type to human-readable string.
     * This is used for generating the usage string and also in 
     * TypeOptionParser to encode the argument type in the options dictionary.
     * @see help 
     * @param type Type.
     * @return Converted string.
     */
    std::string type2string(type_t type);
    /**
     * Type of "options" dictionary.
     */
    typedef std::map<std::string, std::string> options_type;
    options_type options;                   //!< Dictionary of options. Use options[key] to access.
    std::vector<std::string> arguments;     //!< List of positional arguments.
protected:
    /**
     * Set option in dictionary.
     * Function the parser uses to store an option argument or a default value in 
     * options dictionary. Can be overridden in a subclass to store additional info.
     * @param option Option object.
     * @param argument Option argument to store.
     */
    virtual void set_option(Option& option, std::string argument);

private:
    std::string use_msg;            //!< Usage message.
    std::vector<Option> opts;       //!< List of options the parser knows about.

    /* Helper for finding which Option (if any) argv[index] matches. */
    void find_opt_short(int argc, char **argv, int &index);
    void find_opt_long(int argc, char **argv, int &index);
};

 template<typename T> T OptionParser::get(std::string option_name)
   {
     std::string s = get_option(option_name);
     std::istringstream strin(s);
     T t;
     strin >> t;
     return t;
   }

} // namespace optparse

#endif
