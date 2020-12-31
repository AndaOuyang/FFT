#include "fft.hpp"

// nooooo, there is no need for complex_number class, because there is already a std::complex<value_type>





namespace fft{
        // class complext number
    class complex_number {
    public:
        // --------------------------- rounding error threshold ------------------------------------
        static double constexpr epsilon = 0.0000001;

        
        //----------------------------constructors---------------------------------
        /*
        * default constructor, return 0+j0
        * 
        */
        complex_number() noexcept;

        /*
        * construct by the real value
        * 
        * @param real: the real component of this complex number 
        */
        complex_number(double real) noexcept;

        /*
        * construct by the real value and imaginary value
        * 
        * @param real: the real component of this complex number 
        * @param imag: the imaginary component of this complex number
        */
        complex_number(double real, double imag);

        /*
        * copy constructor
        * 
        */
        complex_number(complex_number const& oth) noexcept;

        /*
        * fake move constructor. complex_number only has 2 double attributs, not memory intensive. just copy the 2 component
        * 
        */
        complex_number(complex_number&& oth) noexcept;

        //---------------------------destructor------------------------------------
        ~complex_number() = default;

        //---------------------------operators-------------------------------------
        auto operator=(complex_number const&) noexcept -> complex_number&;
        auto operator=(complex_number&&) noexcept -> complex_number&;
        auto operator+() const noexcept -> complex_number;
        auto operator-() const noexcept -> complex_number;
        auto operator+=(complex_number const&) -> complex_number&;
        auto operator-=(complex_number const&) -> complex_number&;
        auto operator*=(double) noexcept -> complex_number&;
        auto operator/=(double) -> complex_number&;
        auto operator*=(complex_number const&) noexcept -> complex_number&;
        auto operator/=(complex_number const&) -> complex_number&;


        //--------------------------friends----------------------------------------
		friend auto operator==(complex_number const&, complex_number const&) noexcept -> bool;
		friend auto operator!=(complex_number const&, complex_number const&) noexcept -> bool;
		friend auto operator+(complex_number const& lhs, complex_number const& rhs)
		   -> complex_number;
		friend auto operator-(complex_number const& lhs, complex_number const& rhs)
		   -> complex_number;
		friend auto operator*(complex_number const&, double) noexcept -> complex_number;
		friend auto operator/(complex_number const&, double) -> complex_number;
		friend auto operator<<(std::ostream&, complex_number const&) noexcept -> std::ostream&;

		//----------------------Utility functions----------------------------------
        // non const
        friend auto real(complex_number& v) -> double;
        friend auto imag(complex_number& v) -> double;
        
        // const
        friend auto real(complex_number const& v) -> double;
        friend auto imag(complex_number const& v) -> double;
        friend auto argument(complex_number const& v) -> double;
        friend auto phase(complex_number const& v) -> double;
        friend auto abs(complex_number const& v) -> double;

        //-----------------------member functions----------------------------------
        // non const
        [[nodiscard]] auto real() -> double&;
        [[nodiscard]] auto imag() -> double&;

        // const
        [[nodiscard]] auto real() const -> double;
        [[nodiscard]] auto imag() const -> double;
        [[nodiscard]] auto argument() const -> double;
        [[nodiscard]] auto phase() const -> double;
		[[nodiscard]] auto abs() const -> double;
    
    private:
        // ---------------------------------- artributes -------------------------------------------
        double real_;
        double imag_;
    }; // class complex number

    // --------------------------- utility functions of complex_number -----------------------------
    auto real(complex_number const& v) -> double;
    auto imag(complex_number const& v) -> double;
    auto argument(complex_number const& v) -> double;
    auto phase(complex_number const& v) -> double;
    auto abs(complex_number const& v) -> double;

    // ----------------------------------- fft -----------------------------------------------------


    // class complext number
    //----------------------------constructors---------------------------------
    complex_number::complex_number() noexcept: real_{0}, imag_{0} {}
    complex_number::complex_number(double real): real_{real}, imag_{0} {}
    complex_number::complex_number(double real, double imag): real_{real}, imag_{0} {}
    complex_number::complex_number(complex_number const& oth) noexcept: real_{oth.real_}, imag_{oth.imag_} {}
    complex_number::complex_number(complex_number&& oth) noexcept: real_{oth.real_}, imag_{oth.imag_} {}

    //---------------------------operators-------------------------------------
    auto complex_number::operator=(complex_number const& oth) noexcept -> complex_number&{
        if (this == &oth){
            return *this;
        }
        this->real_ = oth.real_;
        this->imag_ = oth.imag_;
        return *this;
    }
    auto complex_number::operator=(complex_number&& oth) noexcept -> complex_number&{
        // do not have much meomory to move, why not just copy the 2 doubles?
        *this = oth;
        return *this;
    }
    auto complex_number::operator+() const noexcept -> complex_number{
        return *this;
    }
    auto complex_number::operator-() const noexcept -> complex_number{
        return complex_number(-this->real_, -this->imag_);
    }
    auto complex_number::operator+=(complex_number const& oth) -> complex_number&{
        this->real_ += oth.real_;
        this->imag_ += oth.imag_;
        return *this;
    }
    auto complex_number::operator-=(complex_number const& oth) -> complex_number&{
        this->real_ -= oth.real_;
        this->imag_ -= oth.imag_;
        return *this;
    }
    auto complex_number::operator*=(double factor) noexcept -> complex_number&{
        this->real_ *= factor;
        this->imag_ *= factor;
        return *this;
    }
    auto complex_number::operator/=(double dividend) -> complex_number&{
        this->real_ /= dividend;
        this->imag_ /= dividend;
        return *this;
    }
    auto complex_number::operator*=(complex_number const& oth) noexcept -> complex_number&{
        this->real_ = this->real_ * oth.real_ - this->imag_ * oth.imag_;
        this->imag_ = this->real_ * oth.imag_ + this->imag_ * oth.real_;
        return *this;
    }
    auto complex_number::operator/=(complex_number const& oth) -> complex_number&{

    }
}
