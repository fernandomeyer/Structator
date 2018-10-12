################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/afAlphabet.c \
../src/afChaining.c \
../src/afComputeAflk.c \
../src/afCons.c \
../src/afConsAfLink.c \
../src/afResultsProcessing.c \
../src/afSearch.c \
../src/afSearchStr.c \
../src/afStreamHandling.c \
../src/construct.c 

OBJS += \
./src/afAlphabet.o \
./src/afChaining.o \
./src/afComputeAflk.o \
./src/afCons.o \
./src/afConsAfLink.o \
./src/afResultsProcessing.o \
./src/afSearch.o \
./src/afSearchStr.o \
./src/afStreamHandling.o \
./src/construct.o 

C_DEPS += \
./src/afAlphabet.d \
./src/afChaining.d \
./src/afComputeAflk.d \
./src/afCons.d \
./src/afConsAfLink.d \
./src/afResultsProcessing.d \
./src/afSearch.d \
./src/afSearchStr.d \
./src/afStreamHandling.d \
./src/construct.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


