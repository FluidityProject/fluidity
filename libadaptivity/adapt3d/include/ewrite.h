#ifndef EWRITE_H
#define EWRITE_H

#ifndef ewrite
#define ewrite(priority, format) if (priority<=cdl) write(dunit(priority), format) 
#endif

#ifndef EWRITE
#define EWRITE(priority, format) ewrite(priority, format)
#endif

#endif
