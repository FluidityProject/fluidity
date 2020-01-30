/*
 * Implementation of optparse.h
 * 
 * Roman Geus, 2005
 *
 */

#include <iostream>
#include <sstream>
#include "optparse.hh"

using namespace std;

namespace optparse {

/********************** Option **********************/

Option::Option (string shrt, string lng, string dest,
                string hlp, action_t act, string dfault, type_t type, string allowed_args)
        : action (act),
        shrt_flag (shrt),
        lng_flag (lng),
        help (hlp),
        destination (dest),
        dfault_(dfault),
        type_(type),
        allowed_args_(allowed_args)
{}

Option::~Option () {}

bool 
Option::is_allowed(std::string argument) {
    // Test type
    // FIXME: implement that
    
    // Test allowed values
    if (allowed_args_ == "")
        return true;
    
    size_t pos_begin = 0;
    do {
        size_t pos_comma = allowed_args_.find(',', pos_begin);
        if (pos_comma == string::npos)
            pos_comma = allowed_args_.size();
        if (allowed_args_.substr(pos_begin, pos_comma-pos_begin) == argument)
            return true;
        pos_begin = pos_comma + 1;
    } while(pos_begin < allowed_args_.size());    
    return false;
}

/******************** OptionParser *******************/

OptionParser::OptionParser (string usage)
    : use_msg (usage) 
{ }

OptionParser::~OptionParser () {}

void
OptionParser::add_option(string shrt_flag, string lng_flag, string destination,
                         string help, action_t act, type_t type, string dfault,
                         string allowed_values) 
{
    Option option(shrt_flag, lng_flag, destination, help, act, dfault, type, allowed_values);
    
    /* Add the option to our list of options. */
    opts.push_back(option);

    /* Set the default value for this option, this not only allows you to
     * set default values but insures that every option's destination will
     * be in our dictionary, even if the value is only "".
     */
    set_option(option, dfault);
}

string OptionParser::get_option(std::string option_name) {
    string arg(options[option_name]);
    size_t colon_pos = arg.find(':');
    return arg.substr(0, colon_pos);
}

void
OptionParser::parse_args (int argc, char **argv) {
    /* Walk through the arguments and sort them into options
     * or arguments.
     */
    for (int i = 1; i < argc; i++) {
        /* If we start with a '-' we're probably a <flag><value> pair
         * so we need to figure out which option it is. find_opt() is
         * where the real work is done.
         */
        if (argv[i][0] == '-')
            if (argv[i][1] == '-')
                find_opt_long(argc, argv, i);
            else
                find_opt_short(argc, argv, i);

        /* If we're looking at an argument (i.e. a value with no flag who's
         * meaning is determined by position) just append it into the
         * arguments list.
         */
        else
            arguments.insert(arguments.end(), argv[i]);
    }
}

void
OptionParser::help (ostream& os) {
    const size_t WORD_WRAP = 50;      // max width of last column
    
    // Determine column width for short options
    size_t shrt_flag_max_len = 0;
    for (size_t i=0; i < opts.size(); ++ i) {
        stringstream buf;
        buf << opts[i].shrt_flag;
        if (opts[i].shrt_flag != "" && opts[i].action == STORE)
            buf << " " << type2string(opts[i].type_);
        if (buf.str().size() > shrt_flag_max_len)
            shrt_flag_max_len = buf.str().size();
    }
    
    // Determine column width for long options
    size_t lng_flag_max_len = 0;
    for (size_t i=0; i < opts.size(); ++ i) {
        stringstream buf;
        buf << opts[i].lng_flag;
        if (opts[i].action == STORE)
            buf << "=" << type2string(opts[i].type_);
        if (buf.str().size() > lng_flag_max_len)
            lng_flag_max_len = buf.str().size();
    }
    
    os << use_msg << endl;
    for (size_t i=0; i < opts.size(); ++ i) {
        stringstream line;
        line << "  ";
        
        // short option column
        stringstream shrt_buf;
        shrt_buf << opts[i].shrt_flag;
        if (opts[i].shrt_flag != "" && opts[i].action == STORE)
            shrt_buf << " " << type2string(opts[i].type_);
        line << ' ' << shrt_buf.str();
        for (size_t k = 0; k < shrt_flag_max_len-shrt_buf.str().size()+2; ++ k)
            line << ' ';
            
        // long option column
        stringstream buf;
        buf << opts[i].lng_flag;
        if (opts[i].action == STORE)
            buf << "=" << type2string(opts[i].type_);        
        line << buf.str();
        for (size_t k = 0; k < lng_flag_max_len-buf.str().size()+4; ++ k)
            line << ' ';
        
        // help column
        os << line.str();
        size_t help_col = line.str().size();
        line.str("");
        
        line << opts[i].help;
        bool is_allowed = opts[i].allowed_args_ != "";
        bool is_default = opts[i].action == STORE && opts[i].dfault_ != "";
        if (is_allowed || is_default) {
            line << " (";
            if (is_allowed) {
                line << "possible values=";
                size_t pos_begin = 0;
                size_t pos_comma = opts[i].allowed_args_.find(',', pos_begin);
                while (pos_comma != string::npos) {
                    line << '\'' << opts[i].allowed_args_.substr(pos_begin, pos_comma-pos_begin) << "\',";
                    pos_begin = pos_comma + 1;
                    pos_comma = opts[i].allowed_args_.find(',', pos_begin);
                }
                line << '\'' << opts[i].allowed_args_.substr(pos_begin) << '\'';
                if (is_default)
                    line << ", ";
            }
            if (is_default)
                line << "default=\'" << opts[i].dfault_ << "\'";
            line << ')';
        }
        
        // split over several lines
        size_t begin_pos = 0;
        size_t last_pos = 0;
        size_t next_pos;
        
        do {        
            next_pos = line.str().find_first_of(" ,", last_pos+1);
            if (next_pos == string::npos)
                next_pos = line.str().size()-1;
            if (last_pos == line.str().size()-1 || (next_pos+1 - begin_pos > WORD_WRAP && last_pos > begin_pos)) {
                if (begin_pos != 0)
                    for (size_t k = 0; k < help_col+2; ++ k)
                        os << ' ';
                os << line.str().substr(begin_pos, last_pos+1 - begin_pos) << endl;
                begin_pos = last_pos+1;
                last_pos = begin_pos;
            } else
                last_pos = next_pos;
        } while (begin_pos != line.str().size());
    }        
 }

void
OptionParser::find_opt_short(int argc, char **argv, int &index) {
    /* Step through our list of known options. */
    for (size_t i = 0; i < opts.size(); i++) {
        /* Uses the overridden == operator for the Options class
         * to compare argv[index] to the flags of each Option.
         */
        if (opts[i].shrt_flag == (string)argv[index]) {
            switch (opts[i].action) {
            case STORE_FALSE:
                set_option(opts[i], "0");
                break;
            case STORE_TRUE:
                set_option(opts[i], "1");
                break;
            case STORE:
                /* Set the value and return if we've found a match. */
                if (index >= argc-1)
                    throw OptionError(argv[index], "Missing option argument");
                set_option(opts[i], argv[index+1]);
                index++;                
                break;
            default:
                break;
            };
            return;
        }
    }
    
    /* If we haven't found a match this is not a known argument. */
    throw OptionError(argv[index], "Unknown option");
}

void
OptionParser::find_opt_long(int argc, char **argv, int &index) {
    // Split option and argument at "="
    string option;
    string argument;
    string arg_str(argv[index]);
    size_t equal_pos = arg_str.find('=');
    if (equal_pos != string::npos) {
        option = arg_str.substr(0, equal_pos);
        argument = arg_str.substr(equal_pos+1, string::npos);
    } else
        option = arg_str;
    
    /* Step through our list of known options. */
    for (size_t i = 0; i < opts.size(); i++) {
        if (opts[i].lng_flag == option) {
            switch (opts[i].action) {
            case STORE_FALSE:
                if (argument != "")
                    throw OptionError(option, "No argument expected for this option");
                set_option(opts[i], "0");
                break;
            case STORE_TRUE:
                if (argument != "")
                    throw OptionError(option, "No argument expected for this option");
                set_option(opts[i], "1");
                break;
            case STORE:
                if (argument == "")
                    throw OptionError(option, "Missing option argument");
                set_option(opts[i], argument);
                break;
            default:
                break;
            };
            return;
        }
    }
    
    /* If we haven't found a match this is not a known argument. */
    throw OptionError(option, "Unknown option");
}

void OptionParser::set_option(Option& option, std::string argument) {
    if (option.is_allowed(argument))
        options[option.destination] = argument;
    else
        throw OptionError("", "Invalid argument for option");
}

string OptionParser::type2string(type_t type) {
    if (type == STRING)
        return string("STRING");
    else if (type == DOUBLE)
        return string("DOUBLE");
    else if (type == INT)
        return string("INT");
    else if (type == BOOL)
        return string("BOOL");
    else
        return "";
}

template<> bool OptionParser::get<bool>(std::string option_name)
{
  return get_option(option_name) == "1";
}

} // namespace optparse
