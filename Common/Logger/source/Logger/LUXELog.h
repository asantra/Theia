#ifndef LUXELOG_H
#define LUXELOG_H
/* Copyright(c) 1998-1999, LUXECE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TClass.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TString.h>

using std::ostream;

// deprecation macro
#if defined(__GNUC__) && (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
# define LUXEROOT_DEPRECATED(func) func  __attribute__ ((deprecated))
#elif defined(_MSC_VER) && _MSC_VER >= 1300
# define LUXEROOT_DEPRECATED(func) __declspec(deprecated) func
#else
# define LUXEROOT_DEPRECATED(func) func
#endif

/**
 * class for logging debug, info and error messages
 */
class LUXELog: public TObject
{
 public:

		// Log4j log levels: TRACE, DEBUG, INFO, WARN, ERROR, FATAL
  enum EType_t {kFatal = 0, kError, kWarning, kInfo, kDebug, kMaxType};
  typedef void (*LUXELogNotification)(EType_t type, const char* message );

		// NB: singleton constructor & destructor should not be public!
		// LUXEROOT_DEPRECATED(LUXELog());
		// LUXEROOT_DEPRECATED(virtual ~LUXELog());

		// NB: singleton deprecated static instance method
		// LUXEROOT_DEPRECATED(static LUXELog* Instance() {return fgInstance;};)

		// get root logger singleton instance
		static LUXELog *GetRootLogger();

		// delete root logger singleton instance
		static void DeleteRootLogger(); 

		// NB: the following functions should not be static
		// NB: deprecated: logging configuration should be made through to a configuration file
  static void  EnableCoreDump(Bool_t enabled);
  static void MakeCoreDump(const char *fout);	
  static void  EnableDebug(Bool_t enabled);
  static void  SetGlobalLogLevel(EType_t type);
  static Int_t GetGlobalLogLevel();
  static void  SetGlobalDebugLevel(Int_t level);
  static Int_t GetGlobalDebugLevel();
  static void  SetModuleDebugLevel(const char* module, Int_t level);
  static void  ClearModuleDebugLevel(const char* module);
  static void  SetClassDebugLevel(const char* className, Int_t level);
  static void  ClearClassDebugLevel(const char* className);

  static void  SetStandardOutput();
  static void  SetStandardOutput(EType_t type);
  static void  SetErrorOutput();
  static void  SetErrorOutput(EType_t type);
  static void  SetFileOutput(const char* fileName);
  static void  SetFileOutput(EType_t type, const char* fileName);
  static void  SetStreamOutput(ostream* stream);
  static void  SetStreamOutput(EType_t type, ostream* stream);
  static void  SetLogNotification(LUXELogNotification pCallBack);
  static void  SetLogNotification(EType_t type, LUXELogNotification pCallBack);
  static void  Flush();

  static void  SetHandleRootMessages(Bool_t on);

  static void  SetPrintType(Bool_t on);
  static void  SetPrintType(EType_t type, Bool_t on);
  static void  SetPrintModule(Bool_t on);
  static void  SetPrintModule(EType_t type, Bool_t on);
  static void  SetPrintScope(Bool_t on);
  static void  SetPrintScope(EType_t type, Bool_t on);
  static void  SetPrintLocation(Bool_t on);
  static void  SetPrintLocation(EType_t type, Bool_t on);

  static void  SetPrintRepetitions(Bool_t on);

  static void  WriteToFile(const char* name, Int_t option = 0);

  // the following public methods are used by the preprocessor macros 
  // and should not be called directly
  static Bool_t IsDebugEnabled() {return fgDebugEnabled;}
  static Int_t GetDebugLevel(const char* module, const char* className);
  static void  Message(UInt_t level, const char* message, 
                       const char* module, const char* className,
                       const char* function, const char* file, Int_t line);
  static void  Debug(UInt_t level, const char* message, 
                     const char* module, const char* className,
                     const char* function, const char* file, Int_t line);

  static Int_t RedirectStdoutTo(EType_t type, UInt_t level, const char* module, 
                                const char* className, const char* function,
                                const char* file, Int_t line, Bool_t print);
  static Int_t RedirectStderrTo(EType_t type, UInt_t level, const char* module, 
                                const char* className, const char* function,
                                const char* file, Int_t line, Bool_t print);
  static void  RestoreStdout(Int_t original);
  static void  RestoreStderr(Int_t original);

  static ostream& Stream(EType_t type, UInt_t level,
                         const char* module, const char* className,
                         const char* function, const char* file, Int_t line);
  static void TestException(Int_t level=10); 
 private:

		// constructor is made private for implementing a singleton
		LUXELog();
		virtual ~LUXELog();

		// NOT IMPLEMENTED?
  LUXELog(const LUXELog& log);
  LUXELog& operator = (const LUXELog& log);

  void           ReadEnvSettings();

  static void    RootErrorHandler(Int_t level, Bool_t abort, 
				  const char* location, const char* message);

  void           CloseFile(Int_t type);
  FILE*          GetOutputStream(Int_t type);

  UInt_t         GetLogLevel(const char* module, const char* className) const;
  void           PrintMessage(UInt_t type, const char* message, 
                              const char* module, const char* className,
                              const char* function, 
                              const char* file, Int_t line);

  void           PrintString(Int_t type, FILE* stream, const char* format, ...);
  void           PrintRepetitions();

  Int_t          RedirectTo(FILE* stream, EType_t type, UInt_t level,
                            const char* module, const char* className,
                            const char* function,
                            const char* file, Int_t line, Bool_t print);

  ostream&       GetStream(EType_t type, UInt_t level,
                           const char* module, const char* className,
                           const char* function, const char* file, Int_t line);

  enum {kDebugOffset = kDebug-1};

  static LUXELog* fgInstance;                 //! pointer to current instance
  static Bool_t  fgDebugEnabled;             // flag for debug en-/disabling
  static Bool_t  fgCoreEnabled;             // flag for core dump en-/disabling

  UInt_t         fGlobalLogLevel;            // global logging level
  TObjArray      fModuleDebugLevels;         // debug levels for modules
  TObjArray      fClassDebugLevels;          // debug levels for classes

  Int_t          fOutputTypes[kMaxType];     // types of output streams
  TString        fFileNames[kMaxType];       // file names
  FILE*          fOutputFiles[kMaxType];     //! log output files
  ostream*       fOutputStreams[kMaxType];   //! log output streams

  Bool_t         fPrintType[kMaxType];       // print type on/off
  Bool_t         fPrintModule[kMaxType];     // print module on/off
  Bool_t         fPrintScope[kMaxType];      // print scope/class name on/off
  Bool_t         fPrintLocation[kMaxType];   // print file and line on/off

  Bool_t         fPrintRepetitions;          // print number of repetitions instead of repeated message on/off

  Int_t          fRepetitions;               //! counter of repetitions
  UInt_t         fLastType;                  //! type of last message
  TString        fLastMessage;               //! last message
  TString        fLastModule;                //! module name of last message
  TString        fLastClassName;             //! class name of last message
  TString        fLastFunction;              //! function name of last message
  TString        fLastFile;                  //! file name of last message
  Int_t          fLastLine;                  //! line number of last message
  LUXELogNotification fCallBacks[kMaxType];   //! external notification callback

  ClassDef(LUXELog, 1)   // class for logging debug, info and error messages
};


// module name macro
#ifdef _MODULE_
# define MODULENAME() _MODULE_
#else
# define MODULENAME() "NoModule"
#endif

// function name macro
#if defined(__GNUC__) || defined(__ICC) || defined(__ECC) || defined(__APPLE__)
# define FUNCTIONNAME() __FUNCTION__
// #elif defined(__HP_aCC) || defined(__alpha) || defined(__DECCXX)
// #define FUNCTIONNAME() __FUNC__
#else
# define FUNCTIONNAME() "???"
#endif

// redirection
/** 
 * Redirect output to std::cout to specified log stream 
 * 
 * @param type      Type of stream to re-direct to
 * @param level     Level of output
 * @param scope     Scope
 * @param whatever  Any code that will output to std::cout 
 */
#define REDIRECTSTDOUT(type, level, scope, whatever) \
  do {Int_t originalStdout = LUXELog::RedirectStdoutTo(type, level, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__, kFALSE); \
    whatever; LUXELog::RestoreStdout(originalStdout);} while(false)
/** 
 * Redirect output to std::cerr to specified log stream 
 * 
 * @param type      Type of stream to re-direct to
 * @param level     Level of output
 * @param scope     Scope
 * @param whatever  Any code that will output to std::cout 
 */
#define REDIRECTSTDERR(type, level, scope, whatever) \
  do {Int_t originalStderr = LUXELog::RedirectStderrTo(type, level, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__, kFALSE); \
    whatever; LUXELog::RestoreStderr(originalStderr);} while(false)
/** 
 * Redirect output to std::cout and std::cerr to specified log stream 
 * 
 * @param type      Type of stream to re-direct to
 * @param level     Level of output
 * @param scope     Scope
 * @param whatever  Any code that will output to std::cout or std::cerr
 */
#define REDIRECTSTDOUTANDSTDERR(type, level, scope, whatever) \
  do {Int_t originalStdout = LUXELog::RedirectStdoutTo(type, level, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__, kFALSE); \
    Int_t originalStderr = LUXELog::RedirectStderrTo(type, level, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__, kFALSE); \
    whatever; LUXELog::RestoreStderr(originalStderr); LUXELog::RestoreStdout(originalStdout);} while(false)


// debug level
#ifdef LOG_NO_DEBUG
# define LUXEDebugLevel() -1
# define LUXEDebugLevelClass() -1
# define LUXEDebugLevelGeneral(scope) -1
#else
/** 
 * Get the object scope debug level
 */
# define LUXEDebugLevel() ((LUXELog::IsDebugEnabled()) ? LUXELog::GetDebugLevel(MODULENAME(), ClassName()) : -1)
/** 
 * Get the class (ROOT-enabled) debug level
 */
# define LUXEDebugLevelClass() ((LUXELog::IsDebugEnabled()) ? LUXELog::GetDebugLevel(MODULENAME(), Class()->GetName()) : -1)
/**
 * Get the debug level associated with scope 
 * @param scope Scope 
 */
# define LUXEDebugLevelGeneral(scope) ((LUXELog::IsDebugEnabled()) ? LUXELog::GetDebugLevel(MODULENAME(), scope) : -1)
#endif

// debug messages
#ifdef LOG_NO_DEBUG
# define LUXEDebug(level, message) do { } while (false)
# define LUXEDebugClass(level, message) do { } while (false)
# define LUXEDebugGeneral(scope, level, message) do { } while (false)
# define LUXEDebugF(level, message,...) do { } while (false)
# define LUXEDebugClassF(level, message,...) do { } while (false)
# define LUXEDebugGeneralF(scope, level, message,...) do { } while (false)
#else

// inspired by log4cxx code (see log4cxx/Logger.h)
// implements GCC branch prediction for increasing logging performance
# if !defined(LUXEROOT_UNLIKELY)
#  if defined(__GNUC__) && (__GNUC__ >= 3)
/**
 * Provides optimization hint to the compiler
 * to optimize for the expression being false.
 * @param expr boolean expression.
 * @returns value of expression.
 */
#   define LUXEROOT_UNLIKELY(expr) __builtin_expect(expr, 0)
#  else
/**
 * Provides optimization hint to the compiler
 * to optimize for the expression being false.
 * @param expr boolean expression.
 * @returns value of expression.
 */
#   define LUXEROOT_UNLIKELY(expr) expr
#  endif
# endif 

/**
 * 
 * Logs a message to a specified logger with the DEBUG level.
 * 
 * @param logLevel the debug level.
 * @param message message to print in the following format: Form(message).
 * Note, that message should contain balanced parenthesis, like 
 * <code>LUXEDebug(1, Form("Failed to decode line %d of %s", line, filename));</code> 
 */
# define LUXEDebug(logLevel, message) \
        do { if (LUXEROOT_UNLIKELY(LUXELog::IsDebugEnabled() && LUXELog::GetDebugLevel(MODULENAME(), ClassName()) >= logLevel)) {\
	  LUXELog::Debug(logLevel, message, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__); }} while (0)
/**
 * 
 * Logs a message to a specified logger with the DEBUG level.  For use
 * in static member functions of a class 
 * 
 * @param logLevel the debug level.
 * @param message message to print in the following format: Form(message).
 * Note, that message should contain balanced parenthesis, like 
 * <code>LUXEDebug(1, Form("Failed to decode line %d of %s", line, filename));</code> 
 */
# define LUXEDebugClass(logLevel, message) \
	do { if (LUXEROOT_UNLIKELY(LUXELog::IsDebugEnabled() && LUXELog::GetDebugLevel(MODULENAME(), Class()->GetName()) >= logLevel)) {\
	  LUXELog::Debug(logLevel, message, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__); }} while (0)

/**
 * Logs a message to a specified logger with the DEBUG level.  For use
 * in non-ROOT-enabled-class scope.
 * 
 * @param scope the logging scope.
 * @param logLevel the debug level.
 * @param message message to print in the following format: Form(message).
 * Note, that message should contain balanced parenthesis, like 
 * <code>LUXEDebug(1, Form("Failed to decode line %d of %s", line, filename));</code> 
*/
# define LUXEDebugGeneral(scope, logLevel, message) \
	do { if (LUXEROOT_UNLIKELY(LUXELog::IsDebugEnabled() && LUXELog::GetDebugLevel(MODULENAME(), scope) >= logLevel)) {\
	  LUXELog::Debug(logLevel, message, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__); }} while (0)
/** 
 * Macro to output debugging information.  This excepts a printf-like
 * format statement.   Note, at least 3 arguments (in total) must be
 * passed. 
 * 
 * @param logLevel Debug level
 * @param format   Printf-like format. 
 */
# define LUXEDebugF(logLevel,format,...) \
do { if (LUXEROOT_UNLIKELY(LUXELog::IsDebugEnabled() && LUXELog::GetDebugLevel(MODULENAME(), ClassName()) >= logLevel)) { \
    TString m;m.Form(format,__VA_ARGS__);					\
    LUXELog::Debug(logLevel, m, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__); }} while (0)
/** 
 * Outut debug information, filtered on debug level.  For use in
 * static member function of a ROOT-enabled class. This excepts a
 * printf-like format statement.  Note, at least 3 arguments (in
 * total) must be passed.
 * 
 * @param logLevel Debug level
 * @param format   Printf-like format 
 * 
 * @return 
 */
# define LUXEDebugClassF(logLevel,format,...) \
  do { if (LUXEROOT_UNLIKELY(LUXELog::IsDebugEnabled() && LUXELog::GetDebugLevel(MODULENAME(), Class()->GetName()) >= logLevel)) { \
      TString m;m.Form(format,__VA_ARGS__);					\
      LUXELog::Debug(logLevel, m, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__); }} while (0)
/** 
 * Outut debug information, filtered on debug level.  For use in
 * static member function of a non-ROOT-enabled class-scope. This
 * excepts a printf-like format statement.  Note, at least 3 arguments
 * (in total) must be passed.
 * 
 * @param scope    Scope 
 * @param logLevel Debug level
 * @param format   Printf-like format 
 * 
 * @return 
 */
# define LUXEDebugGeneralF(scope,logLevel,format,...) \
  do { if (LUXEROOT_UNLIKELY(LUXELog::IsDebugEnabled() && LUXELog::GetDebugLevel(MODULENAME(), scope) >= logLevel)) { \
      TString m;m.Form(format,__VA_ARGS__);					\
      LUXELog::Debug(logLevel, m, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__); }} while (0)
    
#endif

// redirection to debug
#define StdoutToLUXEDebug(level, whatever) REDIRECTSTDOUT(LUXELog::kDebug, level, ClassName(), whatever)
#define StderrToLUXEDebug(level, whatever) REDIRECTSTDERR(LUXELog::kDebug, level, ClassName(), whatever)
#define ToLUXEDebug(level, whatever) REDIRECTSTDOUTANDSTDERR(LUXELog::kDebug, level, ClassName(), whatever)
#define StdoutToLUXEDebugClass(level, whatever) REDIRECTSTDOUT(LUXELog::kDebug, level, Class()->GetName(), whatever)
#define StderrToLUXEDebugClass(level, whatever) REDIRECTSTDERR(LUXELog::kDebug, level, Class()->GetName(), whatever)
#define ToLUXEDebugClass(level, whatever) REDIRECTSTDOUTANDSTDERR(LUXELog::kDebug, level, Class()->GetName(), whatever)
#define StdoutToLUXEDebugGeneral(scope, level, whatever) REDIRECTSTDOUT(LUXELog::kDebug, level, scope, whatever)
#define StderrToLUXEDebugGeneral(scope, level, whatever) REDIRECTSTDERR(LUXELog::kDebug, level, scope, whatever)
#define ToLUXEDebugGeneral(scope, level, whatever) REDIRECTSTDOUTANDSTDERR(LUXELog::kDebug, level, scope, whatever)

// debug stream objects
#define LUXEDebugStream(level) LUXELog::Stream(LUXELog::kDebug, level, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define LUXEDebugClassStream(level) LUXELog::Stream(LUXELog::kDebug, level, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define LUXEDebugGeneralStream(scope, level) LUXELog::Stream(LUXELog::kDebug, level, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__)


/** 
 * Macro that will output stuff using the logging facilities. 
 * 
 * @param lvl     Message level 
 * @param message Message to show 
 */
#define LUXEMessage(lvl,message) do { \
      LUXELog::Message(lvl, message, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__);} while(false) 
/** 
 * Macro that will output stuff using the logging facilities. 
 * For use in static member function of ROOT-enabled class-scope.
 *
 * @param lvl     Message level 
 * @param message Message to show 
 */
#define LUXEMessageClass(lvl,message) do { \
    LUXELog::Message(lvl, message, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__);} while(false) 
/** 
 * Macro that will output stuff using the logging facilities. 
 * For use in non-ROOT-enabled class-scope.
 *
 * @param scope   Scope 
 * @param lvl     Message level 
 * @param message Message to show 
 */
#define LUXEMessageGeneral(scope,lvl,message) do {			\
    LUXELog::Message(lvl, message, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__);} while(false) 
/** 
 * Print a message using the LUXELog logging facility. This macro
 * accepts printf-like format arguments.  Note, at least 3 arguments
 * must be passed.  
 * @code
 *   LUXEMessageF(1, "foo");        // <-- Failes
 *   LUXEMessageF(1, "foo %d", 42); // <-- OK
 * @endcode
 *
 * @param lvl     Message level
 * @param format  printf-like format
 */
#define LUXEMessageF(lvl,format,...) do { \
  TString m; m.Form(format,__VA_ARGS__); \
  LUXELog::Message(lvl, m, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__);} while(false) 
/** 
 * Print a message using the LUXELog logging facility. This macro
 * accepts printf-like format arguments.  Note, at least 3 arguments
 * must be passed.  
 * @code
 *   LUXEMessageF(1, "foo");        // <-- Failes
 *   LUXEMessageF(1, "foo %d", 42); // <-- OK
 * @endcode
 *
 * This is for static member function in for ROOT-enabled class-scope
 *
 * @param lvl     Message level
 * @param format  printf-like format
 */
#define LUXEMessageClassF(lvl,format,...) do { \
  TString m; m.Form(format,__VA_ARGS__); \
  LUXELog::Message(lvl, m, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__);} while(false) 
/** 
 * Print a message using the LUXELog logging facility. This macro
 * accepts printf-like format arguments.  Note, at least 3 arguments
 * must be passed.  
 * @code
 *   LUXEMessageF(1, "foo");        // <-- Failes
 *   LUXEMessageF(1, "foo %d", 42); // <-- OK
 * @endcode
 *
 * This is for non-ROOT-enabled class-scope
 *
 * @param scope   Scope 
 * @param lvl     Message level
 * @param format  printf-like format
 */
#define LUXEMessageGeneralF(scope,lvl,format,...) do {	\
  TString m; m.Form(format,__VA_ARGS__); \
  LUXELog::Message(lvl, m, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__);} while(false) 

// info messages
#ifdef LOG_NO_INFO
# define LUXEInfo(message) do { } while (false)
# define LUXEInfoClass(message) do { } while (false)
# define LUXEInfoGeneral(scope, message) do { } while (false)
# define LUXEInfoF(message,...) do { } while (false)
# define LUXEInfoClassF(message,...) do { } while (false)
# define LUXEInfoGeneralF(scope, message,...) do { } while (false)
#else
/**
 * Forwards to LUXEMessage with log level of LUXELog::kInfo
 * @see LUXEMessage 
 */
# define LUXEInfo(message)               LUXEMessage(LUXELog::kInfo, message)
/**
 * Forwards to LUXEMessageClass with log level of LUXELog::kInfo
 * @see LUXEMessageClass 
 */
# define LUXEInfoClass(message)          LUXEMessageClass(LUXELog::kInfo, message)
/**
 * Forwards to LUXEMessageGeneral with log level of LUXELog::kInfo
 * @see LUXEMessageGeneral
 */
# define LUXEInfoGeneral(scope, message) LUXEMessageGeneral(scope, LUXELog::kInfo, message)
/**
 * Forwards to LUXEMessageF with log level of LUXELog::kInfo
 * @see LUXEMessageF 
 */
# define LUXEInfoF(message,...)               LUXEMessageF(LUXELog::kInfo, message, __VA_ARGS__)
/**
 * Forwards to LUXEMessageClassF with log level of LUXELog::kInfo
 * @see LUXEMessageClassF 
 */
# define LUXEInfoClassF(message,...)          LUXEMessageClassF(LUXELog::kInfo, message, __VA_ARGS__)
/**
 * Forwards to LUXEMessageGeneralF with log level of LUXELog::kInfo
 * @see LUXEMessageGeneralF
 */
# define LUXEInfoGeneralF(scope,message,...)  LUXEMessageGeneralF(scope, LUXELog::kInfo, message, __VA_ARGS__)
#endif

// redirection to info
#define StdoutToLUXEInfo(whatever) REDIRECTSTDOUT(LUXELog::kInfo, 0, ClassName(), whatever)
#define StderrToLUXEInfo(whatever) REDIRECTSTDERR(LUXELog::kInfo, 0, ClassName(), whatever)
#define ToLUXEInfo(whatever) REDIRECTSTDOUTANDSTDERR(LUXELog::kInfo, 0, ClassName(), whatever)
#define StdoutToLUXEInfoClass(whatever) REDIRECTSTDOUT(LUXELog::kInfo, 0, Class()->GetName(), whatever)
#define StderrToLUXEInfoClass(whatever) REDIRECTSTDERR(LUXELog::kInfo, 0, Class()->GetName(), whatever)
#define ToLUXEInfoClass(whatever) REDIRECTSTDOUTANDSTDERR(LUXELog::kInfo, 0, Class()->GetName(), whatever)
#define StdoutToLUXEInfoGeneral(scope, whatever) REDIRECTSTDOUT(LUXELog::kInfo, 0, scope, whatever)
#define StderrToLUXEInfoGeneral(scope, whatever) REDIRECTSTDERR(LUXELog::kInfo, 0, scope, whatever)
#define ToLUXEInfoGeneral(scope, whatever) REDIRECTSTDOUTANDSTDERR(LUXELog::kInfo, 0, scope, whatever)

// info stream objects
#define LUXEInfoStream() LUXELog::Stream(LUXELog::kInfo, 0, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define LUXEInfoClassStream() LUXELog::Stream(LUXELog::kInfo, 0, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define LUXEInfoGeneralStream(scope) LUXELog::Stream(LUXELog::kInfo, 0, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__)

// warning messages
#ifdef LOG_NO_WARNING
# define LUXEWarning(message) do { } while (false)
# define LUXEWarningClass(message) do { } while (false)
# define LUXEWarningGeneral(scope, message) do { } while (false)
# define LUXEWarningF(message,...) do { } while (false)
# define LUXEWarningClassF(message,...) do { } while (false)
# define LUXEWarningGeneralF(scope, message,...) do { } while (false)
#else
/**
 * Forwards to LUXEMessage with log level of LUXELog::kWarning
 * @see LUXEMessage 
 */
# define LUXEWarning(message)               LUXEMessage(LUXELog::kWarning, message)
/**
 * Forwards to LUXEMessageClass with log level of LUXELog::kWarning
 * @see LUXEMessageClass 
 */
# define LUXEWarningClass(message)          LUXEMessageClass(LUXELog::kWarning, message)
/**
 * Forwards to LUXEMessageGeneral with log level of LUXELog::kWarning
 * @see LUXEMessageGeneral
 */
# define LUXEWarningGeneral(scope, message) LUXEMessageGeneral(scope, LUXELog::kWarning, message)
/**
 * Forwards to LUXEMessageF with log level of LUXELog::kWarning
 * @see LUXEMessageF 
 */
# define LUXEWarningF(message,...)               LUXEMessageF(LUXELog::kWarning, message, __VA_ARGS__)
/**
 * Forwards to LUXEMessageClassF with log level of LUXELog::kWarning
 * @see LUXEMessageClassF 
 */
# define LUXEWarningClassF(message,...)          LUXEMessageClassF(LUXELog::kWarning, message, __VA_ARGS__)
/**
 * Forwards to LUXEMessageGeneralF with log level of LUXELog::kWarning
 * @see LUXEMessageGeneralF
 */
# define LUXEWarningGeneralF(scope,message,...)  LUXEMessageGeneralF(scope, LUXELog::kWarning, message, __VA_ARGS__)
#endif

// redirection to warning
#define StdoutToLUXEWarning(whatever) REDIRECTSTDOUT(LUXELog::kWarning, 0, ClassName(), whatever)
#define StderrToLUXEWarning(whatever) REDIRECTSTDERR(LUXELog::kWarning, 0, ClassName(), whatever)
#define ToLUXEWarning(whatever) REDIRECTSTDOUTANDSTDERR(LUXELog::kWarning, 0, ClassName(), whatever)
#define StdoutToLUXEWarningClass(whatever) REDIRECTSTDOUT(LUXELog::kWarning, 0, Class()->GetName(), whatever)
#define StderrToLUXEWarningClass(whatever) REDIRECTSTDERR(LUXELog::kWarning, 0, Class()->GetName(), whatever)
#define ToLUXEWarningClass(whatever) REDIRECTSTDOUTANDSTDERR(LUXELog::kWarning, 0, Class()->GetName(), whatever)
#define StdoutToLUXEWarningGeneral(scope, whatever) REDIRECTSTDOUT(LUXELog::kWarning, 0, scope, whatever)
#define StderrToLUXEWarningGeneral(scope, whatever) REDIRECTSTDERR(LUXELog::kWarning, 0, scope, whatever)
#define ToLUXEWarningGeneral(scope, whatever) REDIRECTSTDOUTANDSTDERR(LUXELog::kWarning, 0, scope, whatever)

// warning stream objects
#define LUXEWarningStream() LUXELog::Stream(LUXELog::kWarning, 0, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define LUXEWarningClassStream() LUXELog::Stream(LUXELog::kWarning, 0, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define LUXEWarningGeneralStream(scope) LUXELog::Stream(LUXELog::kWarning, 0, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__)


// error messages
/**
 * Forwards to LUXEMessage with log level of LUXELog::kError
 * @see LUXEMessage 
 */
#define LUXEError(message)               LUXEMessage(LUXELog::kError, message)
/**
 * Forwards to LUXEMessageClass with log level of LUXELog::kError
 * @see LUXEMessageClass 
 */
#define LUXEErrorClass(message)          LUXEMessageClass(LUXELog::kError, message)
/**
 * Forwards to LUXEMessageGeneral with log level of LUXELog::kError
 * @see LUXEMessageGeneral
 */
#define LUXEErrorGeneral(scope, message) LUXEMessageGeneral(scope, LUXELog::kError, message)
/**
 * Forwards to LUXEMessageF with log level of LUXELog::kError
 * @see LUXEMessageF 
 */
#define LUXEErrorF(message,...)               LUXEMessageF(LUXELog::kError, message, __VA_ARGS__)
/**
 * Forwards to LUXEMessageClassF with log level of LUXELog::kError
 * @see LUXEMessageClassF 
 */
#define LUXEErrorClassF(message,...)          LUXEMessageClassF(LUXELog::kError, message, __VA_ARGS__)
/**
 * Forwards to LUXEMessageGeneralF with log level of LUXELog::kError
 * @see LUXEMessageGeneralF
 */
#define LUXEErrorGeneralF(scope,message,...)  LUXEMessageGeneralF(scope, LUXELog::kError, message, __VA_ARGS__)

// redirection to error
#define StdoutToLUXEError(whatever) REDIRECTSTDOUT(LUXELog::kError, 0, ClassName(), whatever)
#define StderrToLUXEError(whatever) REDIRECTSTDERR(LUXELog::kError, 0, ClassName(), whatever)
#define ToLUXEError(whatever) REDIRECTSTDOUTANDSTDERR(LUXELog::kError, 0, ClassName(), whatever)
#define StdoutToLUXEErrorClass(whatever) REDIRECTSTDOUT(LUXELog::kError, 0, Class()->GetName(), whatever)
#define StderrToLUXEErrorClass(whatever) REDIRECTSTDERR(LUXELog::kError, 0, Class()->GetName(), whatever)
#define ToLUXEErrorClass(whatever) REDIRECTSTDOUTANDSTDERR(LUXELog::kError, 0, Class()->GetName(), whatever)
#define StdoutToLUXEErrorGeneral(scope, whatever) REDIRECTSTDOUT(LUXELog::kError, 0, scope, whatever)
#define StderrToLUXEErrorGeneral(scope, whatever) REDIRECTSTDERR(LUXELog::kError, 0, scope, whatever)
#define ToLUXEErrorGeneral(scope, whatever) REDIRECTSTDOUTANDSTDERR(LUXELog::kError, 0, scope, whatever)

// error stream objects
#define LUXEErrorStream() LUXELog::Stream(LUXELog::kError, 0, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define LUXEErrorClassStream() LUXELog::Stream(LUXELog::kError, 0, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define LUXEErrorGeneralStream(scope) LUXELog::Stream(LUXELog::kError, 0, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__)


// fatal messages
/**
 * Forwards to LUXEMessage with log level of LUXELog::kFatal
 * @see LUXEMessage 
 */
#define LUXEFatal(message)               LUXEMessage(LUXELog::kFatal, message)
/**
 * Forwards to LUXEMessageClass with log level of LUXELog::kFatal
 * @see LUXEMessageClass 
 */
#define LUXEFatalClass(message)          LUXEMessageClass(LUXELog::kFatal, message)
/**
 * Forwards to LUXEMessageGeneral with log level of LUXELog::kFatal
 * @see LUXEMessageGeneral
 */
#define LUXEFatalGeneral(scope, message) LUXEMessageGeneral(scope, LUXELog::kFatal, message)
/**
 * Forwards to LUXEMessageF with log level of LUXELog::kFatal
 * @see LUXEMessageF 
 */
#define LUXEFatalF(message,...)               LUXEMessageF(LUXELog::kFatal, message, __VA_ARGS__)
/**
 * Forwards to LUXEMessageClassF with log level of LUXELog::kFatal
 * @see LUXEMessageClassF 
 */
#define LUXEFatalClassF(message,...)          LUXEMessageClassF(LUXELog::kFatal, message, __VA_ARGS__)
/**
 * Forwards to LUXEMessageGeneralF with log level of LUXELog::kFatal
 * @see LUXEMessageGeneralF
 */
#define LUXEFatalGeneralF(scope,message,...)  LUXEMessageGeneralF(scope, LUXELog::kFatal, message, __VA_ARGS__)

#endif
